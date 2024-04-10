"""
    CalculateCovariance(d::Design, eigenvalues::Vector{Float64}, gate_eigenvalues::Vector{Float64})

Generates the circuit eigenvalue covariance matrix given the circuit and gate eigenvalues, using the structure of the circuits and design. Be careful to note that if `d.full_covariance` is not true, this will only produce the diagonal of the covariance matrix.
"""
function CalculateCovariance(
    d::Design,
    eigenvalues::Vector{Float64},
    gate_eigenvalues::Vector{Float64},
)
    # Initialise some variables
    T = length(d.tuple_set)
    # Initialise a trivial mapping to check things are working properly
    trivial_pauli = Pauli(Bool[0], 0)
    trivial_row = SparseVector{Int16, Int32}(1, Int32[], Int16[])
    trivial_track = Vector{Int16}[]
    trivial_mapping = Mapping(trivial_pauli, trivial_pauli, trivial_row, trivial_track)
    # Generate the covariance matrix for the circuit eigenvalues
    M = size(d.matrix, 1)
    tuple_times_factor = sum(d.shot_weights .* d.tuple_times)
    covariance = spzeros(Float64, Int32, M, M)
    reentrant_lock = ReentrantLock()
    for i in 1:T
        tuple_covariance_dict = d.covariance_dict_ensemble[i]
        tuple_experiment_number = d.experiment_numbers[i]
        shot_weight = d.shot_weights[i]
        if i == 1
            tuple_offset = 0
        else
            tuple_offset = sum([length(d.mapping_ensemble[j]) for j in 1:(i - 1)])
        end
        # Generate the relevant term in the covariance matrix for each key in the dictionary
        @threads :static for key in collect(keys(tuple_covariance_dict))
            # If the key lies on the diagonal, calculate the variance, otherwise calculate the covariance, which is a little more complicated
            upper_key = CartesianIndex(tuple_offset + key[1], tuple_offset + key[2])
            lower_key = CartesianIndex(tuple_offset + key[2], tuple_offset + key[1])
            if upper_key == lower_key
                @assert key[1] == key[2]
                @assert tuple_covariance_dict[key][1] == trivial_mapping
                # Calculate the eigenvalue
                eigenvalue = eigenvalues[upper_key[1]]
                # Calculate the variance scaling factor
                eigenvalue_experiments = tuple_covariance_dict[key][2]
                scaling_var =
                    (tuple_times_factor / shot_weight) *
                    (tuple_experiment_number / eigenvalue_experiments)
                # Calculate the variance
                if eigenvalue < 1.0
                    eigenvalue_variance = scaling_var * (1 - eigenvalue^2)
                else
                    eigenvalue_variance =
                        scaling_var * (1 - maximum(eigenvalues[eigenvalues .< 1.0])^2)
                end
                # Set the covariance matrix term
                lock(reentrant_lock) do
                    covariance[upper_key] = eigenvalue_variance
                end
            else
                prod_mapping = tuple_covariance_dict[key][1]
                # Calculate the eigenvalues
                eigenvalue_1 = eigenvalues[upper_key[1]]
                eigenvalue_2 = eigenvalues[upper_key[2]]
                eigenvalue_both =
                    exp(-(prod_mapping.design_row' * (-log.(gate_eigenvalues))))
                # Calculate the covariance scaling factor
                key_1 = CartesianIndex(key[1], key[1])
                key_2 = CartesianIndex(key[2], key[2])
                eigenvalue_1_experiments = tuple_covariance_dict[key_1][2]
                eigenvalue_2_experiments = tuple_covariance_dict[key_2][2]
                eigenvalue_both_experiments = tuple_covariance_dict[key][2]
                scaling_cov =
                    (tuple_times_factor / shot_weight) * (
                        (tuple_experiment_number * eigenvalue_both_experiments) /
                        (eigenvalue_1_experiments * eigenvalue_2_experiments)
                    )
                # Calculate the covariance
                eigenvalue_both_covariance =
                    scaling_cov * (eigenvalue_both - eigenvalue_1 * eigenvalue_2)
                if eigenvalue_both_covariance < 0.0
                    @debug "The negative off-diagonal covariance matrix entry $(round(eigenvalue_both_covariance, sigdigits = 3)) has been set to 0; consider taking more shots."
                    eigenvalue_both_covariance = 0.0
                end
                # Set the covariance matrix entries
                lock(reentrant_lock) do
                    covariance[upper_key] = eigenvalue_both_covariance
                    covariance[lower_key] = eigenvalue_both_covariance
                end
            end
        end
    end
    return covariance::SparseMatrixCSC{Float64, Int32}
end

"""
    SyntheticEigenvalues(d::Design)

Generates the circuit eigenvalues, covariance matrix, and gate eigenvalues using the design matrix. Be careful to note that if `d.full_covariance` is not true, this will only produce the diagonal of the covariance matrix.
"""
function SyntheticEigenvalues(d::Design)
    # Get the gate eigenvalues
    gate_eigenvalues = d.code.gate_eigenvalues
    # Generate the circuit eigenvalues
    eigenvalues = exp.(-(d.matrix * (-log.(gate_eigenvalues))))
    # Generate the covariance matrix
    covariance = CalculateCovariance(d, eigenvalues, gate_eigenvalues)
    return (eigenvalues::Vector{Float64}, covariance::SparseMatrixCSC{Float64, Int32})
end

"""
    MeritData(d::Design)

Generate the covariance matrix of the circuit log-eigenvalue and the gate eigenvalues, which are used to calculate the figure of merit.
"""
function MeritData(d::Design)
    # Ensure that the design comes with the full covariance matrix
    @assert d.full_covariance "To calculate figures of merit, the design must have the full covariance matrix."
    # Estimate the covariance matrix of the eigenvalues
    (eigenvalues, covariance) = SyntheticEigenvalues(d)
    eigenvalues_inv_diag = Diagonal(eigenvalues .^ (-1))
    # Use a first-order Taylor approximation to estimate the covariance matrix of the circuit log-eigenvalues
    covariance_log = convert(
        SparseMatrixCSC{Float64, Int},
        Symmetric(eigenvalues_inv_diag * covariance * eigenvalues_inv_diag),
    )
    return covariance_log::SparseMatrixCSC{Float64, Int}
end

#
function SparseCovarianceInv(
    covariance_log::SparseMatrixCSC{Float64, Int},
    mapping_lengths::Vector{Int},
)
    # Use the block diagonal structure of the covariance matrix to speed up computation of the covariance inverse
    mapping_lower = cumsum([1; mapping_lengths[1:(end - 1)]])
    mapping_upper = cumsum(mapping_lengths)
    M = sum(mapping_lengths)
    covariance_log_inv = spzeros(Float64, Int, M, M)
    for i in eachindex(mapping_lengths)
        mapping_range = mapping_lower[i]:mapping_upper[i]
        mapping_length = mapping_lengths[i]
        # The sparse Cholesky decomposition of A actually computes a decomposition of a permuted matrix PAP'=LL'
        # Hence A=P'LL'P, and our GLS factor is the inverse of P'L, L^(-1)P
        # Note also that A^(-1) = P'L^(-1)'L^(-1)P
        block_chol = cholesky(covariance_log[mapping_range, mapping_range])
        block_L_inv = sparse(inv(LowerTriangular(sparse(block_chol.L))))
        block_perm = sparse(1:mapping_length, block_chol.p, ones(mapping_length))
        covariance_log_inv[mapping_range, mapping_range] =
            block_perm' * block_L_inv' * block_L_inv * block_perm
    end
    return covariance_log_inv::SparseMatrixCSC{Float64, Int}
end

"""
    GLSCovariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
"""
function GLSCovariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the inverse of the covariance matrix of the circuit log-eigenvalues
    mapping_lengths = length.(d.mapping_ensemble)
    gate_eigenvalues = d.code.gate_eigenvalues
    gate_eigenvalues_diag = Diagonal(gate_eigenvalues)
    covariance_log_inv = SparseCovarianceInv(covariance_log, mapping_lengths)
    # Calculate the covariance matrix of the gate eigenvalues using a first-order Taylor approximation
    gls_gate_eigenvalues_cov = Symmetric(
        gate_eigenvalues_diag *
        inv(cholesky(Symmetric(Array(d.matrix' * covariance_log_inv * d.matrix)))) *
        gate_eigenvalues_diag,
    )
    return gls_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    WLSCovariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
"""
function WLSCovariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the WLS estimator
    gate_eigenvalues = d.code.gate_eigenvalues
    gate_eigenvalues_diag = Diagonal(gate_eigenvalues)
    covariance_log_diag_inv = Diagonal(covariance_log)^(-1)
    wls_estimator =
        inv(cholesky(Symmetric(Array(d.matrix' * covariance_log_diag_inv * d.matrix)))) *
        d.matrix' *
        covariance_log_diag_inv
    # Calculate the covariance matrix of the gate eigenvalues using a first-order Taylor approximation
    wls_gate_eigenvalues_cov = Symmetric(
        gate_eigenvalues_diag *
        wls_estimator *
        covariance_log *
        wls_estimator' *
        gate_eigenvalues_diag,
    )
    return wls_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

"""
    OLSCovariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
"""
function OLSCovariance(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the OLS estimator
    gate_eigenvalues = d.code.gate_eigenvalues
    gate_eigenvalues_diag = Diagonal(gate_eigenvalues)
    ols_estimator = inv(bunchkaufman(Symmetric(Array(d.matrix' * d.matrix)))) * d.matrix'
    # Calculate the covariance matrix of the gate eigenvalues using a first-order Taylor approximation
    ols_gate_eigenvalues_cov = Symmetric(
        gate_eigenvalues_diag *
        ols_estimator *
        covariance_log *
        ols_estimator' *
        gate_eigenvalues_diag,
    )
    return ols_gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

# 
function LSCovariance(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int},
    ls_type::Symbol,
)
    if ls_type == :gls
        gate_eigenvalues_cov = GLSCovariance(d, covariance_log)
    elseif ls_type == :wls
        gate_eigenvalues_cov = WLSCovariance(d, covariance_log)
    elseif ls_type == :ols
        gate_eigenvalues_cov = OLSCovariance(d, covariance_log)
    else
        throw(error("The estimator type $(ls_type) must be either :gls, :wls, or :ols."))
    end
    return gate_eigenvalues_cov::Symmetric{Float64, Matrix{Float64}}
end

# 
function NRMSEMoments(cov_eigenvalues::Vector{Float64})
    # Calculate the trace of the gate eigenvalue estimator covariance matrix, and its square
    sigma_tr = sum(cov_eigenvalues)
    sigma_sq_tr = sum(cov_eigenvalues .^ 2)
    # Calculate the expectation and variance of the NRMSE
    N = length(cov_eigenvalues)
    expectation = sqrt(sigma_tr) * (1 - sigma_sq_tr / (4 * sigma_tr^2)) / sqrt(N)
    variance = (sigma_sq_tr / (2 * sigma_tr)) * (1 - sigma_sq_tr / (8 * sigma_tr^2)) / N
    return (expectation::Float64, variance::Float64)
end

# 
function GLSMoments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the GLS moments from the GLS gate eigenvalue estimator covariance matrix
    (gls_expectation, gls_variance) =
        NRMSEMoments(eigvals(GLSCovariance(d, covariance_log)))
    return (gls_expectation::Float64, gls_variance::Float64)
end

# 
function WLSMoments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the WLS moments from the WLS gate eigenvalue estimator covariance matrix
    (wls_expectation, wls_variance) =
        NRMSEMoments(eigvals(WLSCovariance(d, covariance_log)))
    return (wls_expectation::Float64, wls_variance::Float64)
end

# 
function OLSMoments(d::Design, covariance_log::SparseMatrixCSC{Float64, Int})
    # Calculate the OLS moments from the OLS gate eigenvalue estimator covariance matrix
    (ols_expectation, ols_variance) =
        NRMSEMoments(eigvals(OLSCovariance(d, covariance_log)))
    return (ols_expectation::Float64, ols_variance::Float64)
end

# 
function LSMoments(
    d::Design,
    covariance_log::SparseMatrixCSC{Float64, Int},
    ls_type::Symbol,
)
    if ls_type == :gls
        (expectation, variance) = GLSMoments(d, covariance_log)
    elseif ls_type == :wls
        (expectation, variance) = WLSMoments(d, covariance_log)
    elseif ls_type == :ols
        (expectation, variance) = OLSMoments(d, covariance_log)
    else
        throw(error("The estimator type $(ls_type) must be either :gls, :wls, or :ols."))
    end
    return (expectation::Float64, variance::Float64)
end

# 
function GLSMerit(d::Design)
    # Calculate the circuit log-eigenvalue covariance matrix and the gate eigenvalues
    covariance_log = MeritData(d)
    # Calculate the pseudoinverse norm and condition number of the design matrix
    design_matrix_singular_vals = zeros(Float64, d.code.N)
    try
        design_matrix_singular_vals = svd(convert(Matrix{Float64}, Array(d.matrix))).S
    catch
        @warn "The default singular value decomposition algorithm failed; falling back to QR iteration."
        design_matrix_singular_vals =
            svd(
                convert(Matrix{Float64}, Array(d.matrix));
                alg = LinearAlgebra.QRIteration(),
            ).S
    end
    design_matrix_norm = design_matrix_singular_vals[1]
    design_matrix_pinv_norm = 1 / design_matrix_singular_vals[end]
    design_matrix_cond_num = design_matrix_norm * design_matrix_pinv_norm
    # Calculate the GLS merit from the GLS gate eigenvalue estimator covariance matrix
    gls_cov_eigenvalues = eigvals(GLSCovariance(d, covariance_log))
    (gls_expectation, gls_variance) = NRMSEMoments(gls_cov_eigenvalues)
    gls_merit = Merit(
        d.tuple_set,
        d.tuple_set_data,
        d.code.code_param,
        d.code.noise_param,
        length(d.code.total_gates),
        d.code.N,
        :gls,
        gls_expectation,
        gls_variance,
        gls_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
    )
    return gls_merit::Merit
end

# 
function WLSMerit(d::Design)
    # Calculate the circuit log-eigenvalue covariance matrix and the gate eigenvalues
    covariance_log = MeritData(d)
    # Calculate the pseudoinverse norm and condition number of the design matrix
    design_matrix_singular_vals = zeros(Float64, d.code.N)
    try
        design_matrix_singular_vals = svd(convert(Matrix{Float64}, Array(d.matrix))).S
    catch
        @warn "The default singular value decomposition algorithm failed; falling back to QR iteration."
        design_matrix_singular_vals =
            svd(
                convert(Matrix{Float64}, Array(d.matrix));
                alg = LinearAlgebra.QRIteration(),
            ).S
    end
    design_matrix_norm = design_matrix_singular_vals[1]
    design_matrix_pinv_norm = 1 / design_matrix_singular_vals[end]
    design_matrix_cond_num = design_matrix_norm * design_matrix_pinv_norm
    # Calculate the WLS merit from the WLS gate eigenvalue estimator covariance matrix
    wls_cov_eigenvalues = eigvals(WLSCovariance(d, covariance_log))
    (wls_expectation, wls_variance) = NRMSEMoments(wls_cov_eigenvalues)
    wls_merit = Merit(
        d.tuple_set,
        d.tuple_set_data,
        d.code.code_param,
        d.code.noise_param,
        length(d.code.total_gates),
        d.code.N,
        :wls,
        wls_expectation,
        wls_variance,
        wls_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
    )
    return wls_merit::Merit
end

# 
function OLSMerit(d::Design)
    # Calculate the circuit log-eigenvalue covariance matrix and the gate eigenvalues
    covariance_log = MeritData(d)
    # Calculate the pseudoinverse norm and condition number of the design matrix
    design_matrix_singular_vals = zeros(Float64, d.code.N)
    try
        design_matrix_singular_vals = svd(convert(Matrix{Float64}, Array(d.matrix))).S
    catch
        @warn "The default singular value decomposition algorithm failed; falling back to QR iteration."
        design_matrix_singular_vals =
            svd(
                convert(Matrix{Float64}, Array(d.matrix));
                alg = LinearAlgebra.QRIteration(),
            ).S
    end
    design_matrix_norm = design_matrix_singular_vals[1]
    design_matrix_pinv_norm = 1 / design_matrix_singular_vals[end]
    design_matrix_cond_num = design_matrix_norm * design_matrix_pinv_norm
    # Calculate the OLS merit from the OLS gate eigenvalue estimator covariance matrix
    ols_cov_eigenvalues = eigvals(OLSCovariance(d, covariance_log))
    (ols_expectation, ols_variance) = NRMSEMoments(ols_cov_eigenvalues)
    ols_merit = Merit(
        d.tuple_set,
        d.tuple_set_data,
        d.code.code_param,
        d.code.noise_param,
        length(d.code.total_gates),
        d.code.N,
        :ols,
        ols_expectation,
        ols_variance,
        ols_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
    )
    return ols_merit::Merit
end

# 
function LSMerit(d::Design, ls_type::Symbol)
    if ls_type == :gls
        merit = GLSMerit(d)
    elseif ls_type == :wls
        merit = WLSMerit(d)
    elseif ls_type == :ols
        merit = OLSMerit(d)
    else
        throw(error("The estimator type $(ls_type) must be either :gls, :wls, or :ols."))
    end
    return merit::Merit
end

# 
function MeritSet(d::Design)
    # Calculate the circuit log-eigenvalue covariance matrix and the gate eigenvalues
    covariance_log = MeritData(d)
    # Calculate the pseudoinverse norm and condition number of the design matrix
    design_matrix_singular_vals = zeros(Float64, d.code.N)
    try
        design_matrix_singular_vals = svd(convert(Matrix{Float64}, Array(d.matrix))).S
    catch
        @warn "The default singular value decomposition algorithm failed; falling back to QR iteration."
        design_matrix_singular_vals =
            svd(
                convert(Matrix{Float64}, Array(d.matrix));
                alg = LinearAlgebra.QRIteration(),
            ).S
    end
    design_matrix_norm = design_matrix_singular_vals[1]
    design_matrix_pinv_norm = 1 / design_matrix_singular_vals[end]
    design_matrix_cond_num = design_matrix_norm * design_matrix_pinv_norm
    # Calculate the GLS merit from the GLS gate eigenvalue estimator covariance matrix
    gls_cov_eigenvalues = eigvals(GLSCovariance(d, covariance_log))
    (gls_expectation, gls_variance) = NRMSEMoments(gls_cov_eigenvalues)
    gls_merit = Merit(
        d.tuple_set,
        d.tuple_set_data,
        d.code.code_param,
        d.code.noise_param,
        length(d.code.total_gates),
        d.code.N,
        :gls,
        gls_expectation,
        gls_variance,
        gls_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
    )
    # Calculate the WLS merit from the WLS gate eigenvalue estimator covariance matrix
    wls_cov_eigenvalues = eigvals(WLSCovariance(d, covariance_log))
    (wls_expectation, wls_variance) = NRMSEMoments(wls_cov_eigenvalues)
    wls_merit = Merit(
        d.tuple_set,
        d.tuple_set_data,
        d.code.code_param,
        d.code.noise_param,
        length(d.code.total_gates),
        d.code.N,
        :wls,
        wls_expectation,
        wls_variance,
        wls_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
    )
    # Calculate the OLS merit from the OLS gate eigenvalue estimator covariance matrix
    ols_cov_eigenvalues = eigvals(OLSCovariance(d, covariance_log))
    (ols_expectation, ols_variance) = NRMSEMoments(ols_cov_eigenvalues)
    ols_merit = Merit(
        d.tuple_set,
        d.tuple_set_data,
        d.code.code_param,
        d.code.noise_param,
        length(d.code.total_gates),
        d.code.N,
        :ols,
        ols_expectation,
        ols_variance,
        ols_cov_eigenvalues,
        design_matrix_cond_num,
        design_matrix_pinv_norm,
        d.tuple_times,
        d.shot_weights,
        d.experiment_numbers,
        d.experiment_number,
    )
    return (gls_merit::Merit, wls_merit::Merit, ols_merit::Merit)
end

# 
function ImhofIntegrand(u::Float64, x::Float64, norm_cov_eigenvalues::Vector{Float64})
    # Integrand of Imhof's method CDF following Eq. 3.2 of `Computing the Distribution of Quadratic Forms in Normal Variables` by J. P. Imhof (1961)
    # Both theta and rho have been simplified as, with reference to Eq. 1.1 we see that h and delta are 1 and 0, respectively
    theta = sum(atan(eigenvalue * u) for eigenvalue in norm_cov_eigenvalues) / 2 - x * u / 2
    rho = prod((1 + (eigenvalue * u)^2) for eigenvalue in norm_cov_eigenvalues)^(1 / 4)
    integrand = sin(theta) / (u * rho)
    return integrand::Float64
end

# 
function nrmse_pdf(
    cov_eigenvalues::Vector{Float64},
    x_values::Vector{Float64},
    epsilon::Float64 = 1e-5,
)
    # Normalise the gate eigenvalue estimator covariance matrix eigenvalues
    N = length(cov_eigenvalues)
    norm_cov_eigenvalues = cov_eigenvalues / N
    # Calculate the moments 
    (expectation, variance) = NRMSEMoments(cov_eigenvalues)
    # Calculate the normal approximation of the PDF
    nrmse_pdf_normal =
        1 / (sqrt(2 * π) * sqrt(variance)) *
        exp.(-(x_values .- expectation) .^ 2 / (2 * variance))
    nrmse_pdf_normal_max = maximum(nrmse_pdf_normal)
    # Perform Imhof's method to calculate the PDF
    # Imhof's method CDF following Eq. 3.2 of `Computing the Distribution of Quadratic Forms in Normal Variables` by J. P. Imhof (1961)
    function ImhofCDF(x)
        return 0.5 - quadgk(u -> ImhofIntegrand(u, x, norm_cov_eigenvalues), 0, Inf)[1] / π
    end
    nrmse_pdf = zeros(length(x_values))
    for (idx, x) in enumerate(x_values)
        # Only perform the calculation for sufficiently large PDF values
        if nrmse_pdf_normal[idx] / nrmse_pdf_normal_max > epsilon
            # Differentiate the CDF to obtain the PDF
            # Note that Imhof's method describes the square, whereas we seek the square root
            nrmse_pdf[idx] = 2 * x * central_fdm(5, 1)(ImhofCDF, x^2)
        end
    end
    return nrmse_pdf::Vector{Float64}
end
