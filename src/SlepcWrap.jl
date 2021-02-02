"""
    I have decided to remove all exclamation marks `!` at the end of routines "modifying" their arguments
    when the name is the same as in SLEPc API. It was too confusing.
"""
module SlepcWrap
    using PetscWrap
    using MPI

    include("load.jl")
    include("init.jl")
    export SlepcInitialize, SlepcFinalize

    include("eps.jl")
    export  SlepcEPS, CEPS,
            EPSCreate,
            EPSDestroy,
            EPSGetConverged,
            EPSGetTolerances,
            EPSGetEigenpair,
            EPSGetEigenvalue,
            EPSSetFromOptions,
            EPSSetOperators,
            EPSSetUp,
            EPSSolve
end