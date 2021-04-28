"""
    Wrapper to SlepcInitializeNoPointers

# Implementation
I don't know if I am supposed to use PetscInt or not...
"""
function SlepcInitialize(args::Vector{String}, filename::String, help::String)
    args2 = ["julia"; args]
    nargs = Cint(length(args2))
    error = ccall( (:SlepcInitializeNoPointers, libslepc),
            PetscErrorCode,
            (Cint,
            Ptr{Ptr{UInt8}},
            Cstring,
            Cstring),
            nargs, args2, filename, help
    )
    @assert iszero(error)
end

SlepcInitialize(args::Vector{String}) = SlepcInitialize(args, "", "")
SlepcInitialize(args::String) = SlepcInitialize(convert(Vector{String}, split(args)), "", "")

"""
    Initialize SLEPc.

If `cmd_line_args == true`, then command line arguments passed to Julia are used as
arguments for SLEPc (leading to a call to `SlepcInitializeNoPointers`).

Otherwise, if `cmd_line_args == false`, initialize SLEPc without arguments (leading
to a call to `SlepcInitializeNoArguments`).
"""
function SlepcInitialize(cmd_line_args::Bool = true)
    if (cmd_line_args)
        SlepcInitialize(ARGS)
    else
        error = ccall((:SlepcInitializeNoArguments, libslepc), PetscErrorCode, ())
        @assert iszero(error)
    end
end

"""
    Wrapper to SlepcFinalize
"""
function SlepcFinalize()
    error = ccall( (:SlepcFinalize, libslepc), PetscErrorCode, ())
    @assert iszero(error)
end