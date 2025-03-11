# Package Performance

QuantumACES has several performance optimisations that are essential for its practical utility.
We now examine some core optimisations that may be of particular interest through the lens of the functions [`optimise_design`](@ref) and [`simulate_aces`](@ref).
These functions optimise experimental designs and simulate ACES noise characterisation experiments, respectively.

First, we discuss the optimisation of experimental designs with [`optimise_design`](@ref).
Crucial to the tractability of this optimisation is the fact that we can perform the optimisation at small scales, such as the syndrome extraction circuits of distance-3 topological codes, and then transfer the optimised design to the syndrome extraction circuits of the same codes at larger distances, while retaining performance characteristics.
As discussed in [arXiv:2404.06545](https://arxiv.org/abs/2404.06545), an ACES experimental design ``(\mathcal{T},\Gamma)`` is parameterised by the tuple set ``\mathcal{T}`` and the shot weights ``\Gamma``.
The tuple set ``\mathcal{T}`` is discrete and cannot be optimised continuously, and so its optimisation relies on repeated calculation of the figure of merit ``\mathcal{F}``, which is given in \cref{eq:aces-figure-of-merit}.
However, the shot weights ``\Gamma`` are continuous, allowing us to optimise them with gradient descent.
This relies on the fast calculation of the gradient of the figure of merit with respect to the shot weights with analytic expressions for these gradients.
Carefully optimised calculation of these analytic expressions is roughly one and a half orders of magnitude faster than forward-mode automatic differentiation, which is much faster than backward-mode automatic differentiation in this case.

We now remark on some considerations in the linear algebra calculations of the figure of merit and its gradient.
The design matrix ``A`` is sparse, and the circuit (log-)eigenvalue estimator covariance matrix ``\Omega`` (``\Omega^\prime``) is block diagonal with the blocks themselves being sparse.
It is important to appropriately leverage this sparsity in calculations.
Matrix multiplications quickly densify sparse matrices, at which point the matrices must be converted to dense form to avoid substantial performance degradation.
In general, QuantumACES calculates matrix inverses for dense symmetric matrices, and so leverages the Cholesky or Bunch-Kaufman decompositions depending on whether the matrix is positive definite or not, respectively.
It is often necessary to invert ``\Omega`` or ``\Omega^\prime``, and this leverages the sparse Cholesky decomposition for separate blocks of the matrix.
Lastly, QuantumACES calculates the condition number and pseudoinverse norm of ``A`` with sparse methods for eigenvalue finding.

Next, we discuss the simulation of ACES noise characterisation experiments with [`simulate_aces`](@ref).
Experiments are simulated with [Stim](https://github.com/quantumlib/Stim), a fast simulator for stabiliser circuits with Pauli noise, with the function [`simulate_stim_estimate`](@ref).
This function simulates experiments with Stim and immediately estimates the circuit eigenvalues measured by the experiment to avoid storing large numbers of measurement outcomes in memory.
By default, Stim outputs measurement outcomes as bits, but can instead pack measurement outcomes into 8-bit unsigned integers, which in practice speeds up sampling of circuit outcomes by roughly a factor of 8.
These bit-packed results can then be processed quickly with careful bit manipulation to estimate the circuit eigenvalues measured by the experiment.
It is also important to automatically split up Stim simulations that attempt to sample too many shots at once to avoid memory issues.
One step in this simulation process that can surprisingly become a bottleneck is the creation of the circuit string for Stim corresponding to each experiment.
If the string manipulation is not optimised, this can become very slow for large circuits.

Once the Stim simulations have been performed to estimate the circuit eigenvalues, we then estimate the gate eigenvalues and Pauli error probabilities with the function [`estimate_gate_noise`](@ref).
We will discuss this process in the context of GLS, though at large scales it is essential to only construct the diagonal of the circuit eigenvalue estimator covariance matrix ``\Omega`` and hence only perform WLS.
First, we calculate the sparse block diagonal Cholesky factorisation of the circuit log-eigenvalue estimator covariance matrix ``\Omega^\prime=LL^\intercal``, and calculate the inverse of the Cholesky factor ``L^{-1}``.
Then we left-divide the scaled design matrix ``L^{-1}A`` by the scaled circuit log-eigenvalues ``L^{-1}\bm{b}`` to obtain the estimated gate log-eigenvalues ``\hat{\bm{x}}``.
The left-division operator in Julia leverages the sparse structure of the inversion problem and is very fast.
Once we have the gate log-eigenvalues, we can straightforwardly estimate the Pauli error probabilities for each gate in the circuit.

However, these estimates are not guaranteed to be valid probability distributions, so we must project the estimates into the probability simplex.
As discussed in [arXiv:2502.21044](https://arxiv.org/abs/2502.21044), we perform this projection in the Mahalanobis distance with the fast convex solver [SCS.jl](https://github.com/jump-dev/SCS.jl).
We can optimise calculation of the precision matrix such that the only matrix inversion required is the highly-optimised inversion of the circuit log-eigenvalue estimator covariance matrix ``\Omega^\prime``.
This projection performs best, at least at small scales, when single-threaded by setting `ENV["OMP_NUM_THREADS"] = "1"` in `~/.julia/config/startup.jl`, but this causes large performance regressions in other linear algebra routines in the package, so the setting is only advised when benchmarking the time taken to perform ACES noise estimation.
