"""
    Wrapper to SlepcInitializeNoArguments
"""
function SlepcInitialize()
    error = ccall((:SlepcInitializeNoArguments, libslepc), PetscErrorCode, ())
    @assert iszero(error)
end

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
    Wrapper to SlepcFinalize
"""
function SlepcFinalize()
    error = ccall( (:SlepcFinalize, libslepc), PetscErrorCode, ())
    @assert iszero(error)
end