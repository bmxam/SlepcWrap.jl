const CEPS = Ptr{Cvoid}

struct SlepcEPS
    ptr::Ref{CEPS}

    SlepcEPS() = new(Ref{Ptr{Cvoid}}())
end

# allows us to pass SlepcEPS objects directly into CEPS ccall signatures
Base.cconvert(::Type{CEPS}, eps::SlepcEPS) = eps.ptr[]

"""
    Wrapper for EPSCreate
"""
function EPSCreate(comm, eps::SlepcEPS)
    error = ccall((:EPSCreate, libslepc), PetscErrorCode, (MPI.MPI_Comm, Ptr{CEPS}), comm, eps.ptr)
    @assert iszero(error)
end

function EPSCreate()
    eps = SlepcEPS()
    EPSCreate(MPI.COMM_WORLD, eps)
    return eps
end

"""
    Wrapper for EPSSetOperators
"""
function EPSSetOperators(eps::SlepcEPS, A::PetscMat, B::PetscMat)
    error = ccall((:EPSSetOperators, libslepc), PetscErrorCode, (CEPS, CMat, CMat), eps, A, B)
    @assert iszero(error)
end

"""
    Wrapper for EPSSetOperators
"""
function EPSSetOperators(eps::SlepcEPS, A::PetscMat)
    error = ccall((:EPSSetOperators, libslepc), PetscErrorCode, (CEPS, CMat, CMat), eps, A, C_NULL)
    @assert iszero(error)
end

"""
    Wrapper for EPSSetup
"""
function EPSSetUp(eps::SlepcEPS)
    error = ccall((:EPSSetUp, libslepc), PetscErrorCode, (CEPS,), eps)
    @assert iszero(error)
end

"""
    Wrapper for EPSSetFromOptions
"""
function EPSSetFromOptions(eps::SlepcEPS)
    error = ccall((:EPSSetFromOptions, libslepc), PetscErrorCode, (CEPS,), eps)
    @assert iszero(error)
end

"""
    Wrapper for EPSSolve
"""
function EPSSolve(eps::SlepcEPS)
    error = ccall((:EPSSolve, libslepc), PetscErrorCode, (CEPS,), eps)
    @assert iszero(error)
end

"""
    Wrapper for EPSGetConverged
"""
function EPSGetConverged(eps::SlepcEPS)
    nconv = Ref{PetscInt}(0)

    error = ccall((:EPSGetConverged, libslepc), PetscErrorCode, (CEPS, Ref{PetscInt}), eps, nconv)
    @assert iszero(error)

    return nconv[]
end

"""
    Wrapper for EPSGetEigenvalue
    `ieig` must be in [1, EPSGetConverged]; it does't start at `0` as in Slepc API.
"""
function EPSGetEigenvalue(eps::SlepcEPS, ieig)
    eigr = Ref{PetscScalar}(0.)
    eigi = Ref{PetscScalar}(0.)

    error = ccall((:EPSGetEigenvalue, libslepc),
                   PetscErrorCode,
                   (CEPS, PetscInt, Ref{PetscScalar}, Ref{PetscScalar}),
                   eps, PetscInt(ieig - 1), eigr, eigi
    )
    @assert iszero(error)

    return eigr[], eigi[]
end

"""
    Wrapper for EPSGetTolerances
"""
function EPSGetTolerances(eps::SlepcEPS)
    tol = Ref{PetscReal}(0)
    maxits = Ref{PetscInt}(0)

    error = ccall((:EPSGetTolerances, libslepc), PetscErrorCode, (CEPS, Ref{PetscReal}, Ref{PetscInt}), eps, tol, maxits)
    @assert iszero(error)

    return tol[], maxits[]
end

"""
    Wrapper for EPSGetEigenpair
    `ieig` must be in [1, EPSGetConverged]; it does't start at `0` as in Slepc API.
"""
function EPSGetEigenpair(eps::SlepcEPS, ieig, vecr::PetscVec, veci::PetscVec)
    eigr = Ref{PetscScalar}(0.)
    eigi = Ref{PetscScalar}(0.)

    error = ccall((:EPSGetEigenpair, libslepc),
                    PetscErrorCode,
                    (CEPS, PetscInt, Ref{PetscScalar}, Ref{PetscScalar}, CVec, CVec),
                    eps, PetscInt(ieig - 1), eigr, eigi, vecr, veci
    )
    @assert iszero(error)

    return eigr[], eigi[], vecr, veci
end

"""
    Wrapper for EPSGetEigenpair without providing pre-allocated vec.
"""
function EPSGetEigenpair(eps::SlepcEPS, mat::PetscMat, ieig)
    eigr = Ref{PetscScalar}(0.)
    eigi = Ref{PetscScalar}(0.)
    vecr, veci = MatCreateVecs(mat)

    return EPSGetEigenpair(eps, ieig, vecr, veci)
end


"""
    Wrapper for EPSDestroy
"""
function EPSDestroy(eps::SlepcEPS)
    error = ccall((:EPSDestroy, libslepc), PetscErrorCode, (Ptr{CEPS},), eps.ptr)
    @assert iszero(error)
end