module SlepcWrap
using PetscWrap
using MPI
using DelimitedFiles

include("const_arch_ind.jl")
# export all items of some enums
for item in Iterators.flatten((
    instances(EPSWhich),
))
    @eval export $(Symbol(item))
end

include("load.jl")
include("init.jl")
export SlepcInitialize, SlepcFinalize

include("eps.jl")
export SlepcEPS, CEPS,
    create_eps,
    get_eig, get_eigs,
    get_eigenvalue, get_eigenvalues,
    get_eigenpair,
    get_tolerances,
    neigs,
    EPSCreate,
    EPSDestroy,
    EPSGetConverged,
    EPSGetTolerances,
    EPSGetEigenpair,
    EPSGetEigenvalue,
    EPSGetEigenvector,
    EPSGetOperators,
    EPSSetFromOptions,
    EPSSetOperators,
    EPSSetTarget,
    EPSSetUp,
    EPSSetWhichEigenpairs,
    EPSSolve,
    EPSView

include("fancy/eps.jl")
export create_eps,
    eigenvalues2file,
    eigenvectors2file,
    get_eig, get_eigs,
    get_eigenvalue, get_eigenvalues,
    get_eigenpair,
    get_tolerances,
    set_target!,
    neigs
end