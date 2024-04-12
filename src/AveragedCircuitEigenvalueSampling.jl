module AveragedCircuitEigenvalueSampling

# Imports
using PythonCall,
    Base.Threads,
    LinearAlgebra,
    SparseArrays,
    Random,
    Distributions,
    Combinatorics,
    StatsBase,
    GLM,
    LsqFit,
    QuadGK,
    FiniteDifferences,
    DataFrames,
    Accessors,
    StructEquality,
    PrettyTables,
    FileIO,
    JLD2

# Abstract types
abstract type AbstractNoiseParameters end
abstract type AbstractCircuitParameters end
abstract type AbstractCircuit end
export AbstractNoiseParameters, AbstractCircuitParameters, AbstractCircuit

# TODO: Clean up exports
# Note that people can still use the functions directly, if for example they do
# import AveragedCircuitEigenvalueSampling as ACES
# ACES.apply_tuple()

# Struct exports
export
    # kwargs.jl
    OptimOptions,
    # tableau.jl structs 
    Tableau,
    Gate,
    Layer,
    # noise.jl structs
    DepolarisingParameters,
    LognormalParameters,
    # circuit.jl structs
    RotatedPlanarParameters,
    UnrotatedPlanarParameters,
    Circuit,
    Code,
    # tuples.jl structs
    TupleSetData,
    # design.jl structs
    Pauli,
    Mapping,
    Design,
    # merit.jl structs
    Merit,
    # scaling.jl structs
    DepolarisingScalingData,
    LognormalScalingData,
    # simulate.jl structs
    ACESData

# Function exports
export
    # tableau.jl functions
    cx!,
    hadamard!,
    phase!,
    x!,
    z!,
    y!,
    cz!,
    sqrt_zz!,
    sqrt_zz_dag!,
    row_sum!,
    measure!,
    reset!,
    apply!,
    make_layer,
    pad_layer,
    unwrap_circuit,
    get_gates,
    label_circuit,
    index_gates,
    apply_tuple,
    # noise.jl functions
    get_gate_probabilities,
    get_gate_eigenvalues,
    # circuit.jl functions
    update_noise,
    get_layer_times,
    rotated_planar_circuit,
    unrotated_planar_circuit,
    get_circuit,
    # tuples.jl functions
    get_basic_tuple_set,
    get_basic_experiment_numbers,
    get_basic_times_harm_mean,
    get_tuple_set_params,
    get_tuple_set_data,
    get_tuple_set,
    # design.jl functions
    get_pauli_prep_set,
    get_prep_layer,
    get_meas_layer,
    calc_mapping,
    calc_mapping_set,
    calc_consistency_set,
    calc_experiment_set,
    calc_covariance_dict,
    get_experiment_layers,
    generate_design,
    complete_design,
    # merit.jl functions
    calc_covariance,
    calc_eigenvalues_covariance,
    calc_covariance_log,
    sparse_covariance_inv,
    calc_gls_covariance,
    calc_wls_covariance,
    calc_ols_covariance,
    calc_ls_covariance,
    nrmse_moments,
    calc_gls_moments,
    calc_wls_moments,
    calc_ols_moments,
    calc_ls_moments,
    calc_gls_merit,
    calc_wls_merit,
    calc_ols_merit,
    calc_ls_merit,
    calc_merit_set,
    nrmse_pdf_integrand,
    nrmse_pdf,
    # weights.jl functions
    get_shot_weights_factor,
    get_shot_weights_factor_inv,
    get_shot_weights_local_grad,
    get_merit_grad,
    calc_gls_merit_grad,
    gls_optimise_weights,
    calc_wls_merit_grad,
    wls_optimise_weights,
    calc_ols_merit_grad,
    ols_optimise_weights,
    optimise_weights,
    compare_ls_optimise_weights,
    # optimise.jl functions
    optimal_expectation,
    step_repetitions,
    optimise_repetitions,
    sample_zipf,
    random_tuple,
    grow_design,
    prune_design,
    grow_design_excursion,
    prune_design_excursion,
    optimise_tuple_set,
    optimise_design,
    # scaling.jl functions
    calc_depolarising_scaling_data,
    calc_lognormal_scaling_data,
    # simulate.jl functions
    get_stim_circuit_string,
    stim_sample,
    batch_shots,
    estimate_eigenvalues,
    fgls_estimate_gate_eigenvalues,
    gls_estimate_gate_eigenvalues,
    wls_estimate_gate_eigenvalues,
    ols_estimate_gate_eigenvalues,
    estimate_gate_probabilities,
    simulate_aces,
    # utils.jl functions
    project_simplex,
    get_support,
    wht_matrix,
    pretty_print,
    get_pauli_string,
    get_mapping_string,
    # io.jl functions
    enter_folder,
    exit_folder,
    code_filename,
    noise_filename,
    tuples_filename,
    design_filename,
    dep_scaling_filename,
    log_scaling_filename,
    aces_data_filename,
    save_design,
    load_design,
    delete_design,
    save_scaling,
    load_scaling,
    delete_scaling,
    save_aces,
    load_aces,
    delete_aces

# Include files
include("kwargs.jl")
include("tableau.jl")
include("noise.jl")
include("circuit.jl")
include("tuples.jl")
include("design.jl")
include("merit.jl")
include("weights.jl")
include("optimise.jl")
include("scaling.jl")
include("simulate.jl")
include("utils.jl")
include("io.jl")

# IntelliSence for Julia VSCode does not work, but this hacky trick fixes that
# It convinces the LSP that the following files are part of src
# Source: https://discourse.julialang.org/t/lsp-missing-reference-woes/98231/16
@static if false
    include("../test/runtests.jl")
    include("../test/merit_tests.jl")
    include("../test/design_tests.jl")
    include("../scalable_aces/rot_optimise.jl")
    include("../scalable_aces/rot_scaling.jl")
    include("../scalable_aces/rot_simulate.jl")
    include("../scalable_aces/rot_simulate_big.jl")
    include("../scalable_aces/unrot_optimise.jl")
    include("../scalable_aces/unrot_scaling.jl")
    include("../scalable_aces/unrot_simulate.jl")
    include("../scalable_aces/unrot_simulate_big.jl")
    include("../scalable_aces/unrot_simulate_big.jl")
end

end
