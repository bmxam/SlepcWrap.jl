module SlepcWrap
    using PetscWrap
    using MPI

    include("load.jl")
    include("init.jl")
    export SlepcInitialize, SlepcFinalize

    include("eps.jl")
    export  SlepcEPS, CEPS,
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
            EPSGetOperators,
            EPSSetFromOptions,
            EPSSetOperators,
            EPSSetUp,
            EPSSolve,
            EPSView

    include("fancy/eps.jl")
    export  create_eps,
            get_eig, get_eigs,
            get_eigenvalue, get_eigenvalues,
            get_eigenpair,
            get_tolerances,
            neigs
end