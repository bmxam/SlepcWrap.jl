const CEPS = Ptr{Cvoid}

struct SlepcEPS
    ptr::Ref{CEPS}
    comm::MPI.Comm

    SlepcEPS(comm::MPI.Comm) = new(Ref{Ptr{Cvoid}}(), comm)
end

# allows us to pass SlepcEPS objects directly into CEPS ccall signatures
Base.cconvert(::Type{CEPS}, eps::SlepcEPS) = eps.ptr[]

"""
    EPSCreate(comm::MPI.Comm, eps::SlepcEPS)

Wrapper for `EPSCreate`
"""
function EPSCreate(comm::MPI.Comm=MPI.COMM_WORLD)
    eps = SlepcEPS(comm)
    error = ccall((:EPSCreate, libslepc), PetscErrorCode, (MPI.MPI_Comm, Ptr{CEPS}), comm, eps.ptr)
    @assert iszero(error)
    return eps
end

"""
    EPSSetOperators(eps::SlepcEPS, A::PetscMat, B::PetscMat)

Wrapper for EPSSetOperators with two matrices
"""
function EPSSetOperators(eps::SlepcEPS, A::PetscMat, B::PetscMat)
    error = ccall((:EPSSetOperators, libslepc), PetscErrorCode, (CEPS, CMat, CMat), eps, A, B)
    @assert iszero(error)
end

"""
    EPSSetOperators(eps::SlepcEPS, A::PetscMat)

Wrapper for EPSSetOperators with only one matrix
"""
function EPSSetOperators(eps::SlepcEPS, A::PetscMat)
    error = ccall((:EPSSetOperators, libslepc), PetscErrorCode, (CEPS, CMat, CMat), eps, A, C_NULL)
    @assert iszero(error)
end

"""
    EPSSetUp(eps::SlepcEPS)

Wrapper for EPSSetup
"""
function EPSSetUp(eps::SlepcEPS)
    error = ccall((:EPSSetUp, libslepc), PetscErrorCode, (CEPS,), eps)
    @assert iszero(error)
end

"""
    EPSSetFromOptions(eps::SlepcEPS)

Wrapper for `EPSSetFromOptions`
"""
function EPSSetFromOptions(eps::SlepcEPS)
    error = ccall((:EPSSetFromOptions, libslepc), PetscErrorCode, (CEPS,), eps)
    @assert iszero(error)
end

"""
    EPSSetTarget(eps::SlepcEPS, target::PetscScalar)

Wrapper for `EPSSetTarget`
"""
function EPSSetTarget(eps::SlepcEPS, target::PetscScalar)
    error = ccall((:EPSSetTarget, libslepc), PetscErrorCode, (CEPS, PetscScalar), eps, target)
    @assert iszero(error)
end
EPSSetTarget(eps::SlepcEPS, target::Number) = EPSSetTarget(eps, PetscScalar(target))

"""
    EPSSetWhichEigenpairs(eps::SlepcEPS, which::EPSWhich)

Wrapper for `EPSSetWhichEigenpairs`
"""
function EPSSetWhichEigenpairs(eps::SlepcEPS, which::EPSWhich)
    error = ccall((:EPSSetWhichEigenpairs, libslepc), PetscErrorCode, (CEPS, EPSWhich), eps, which)
    @assert iszero(error)
end

"""
    EPSGetOperators(eps::SlepcEPS)

Wrapper for `EPSGetOperators`
"""
function EPSGetOperators(eps::SlepcEPS)
    A = PetscMat(eps.comm)
    B = PetscMat(eps.comm)

    error = ccall((:EPSGetOperators, libslepc), PetscErrorCode, (CEPS, Ptr{CMat}, Ptr{CMat}), eps, A.ptr, B.ptr)
    @assert iszero(error)

    return A, B
end

"""
    EPSSolve(eps::SlepcEPS)

Wrapper for EPSSolve
"""
function EPSSolve(eps::SlepcEPS)
    error = ccall((:EPSSolve, libslepc), PetscErrorCode, (CEPS,), eps)
    @assert iszero(error)
end

"""
    EPSGetConverged(eps::SlepcEPS)

Wrapper for `EPSGetConverged` : return the number of converged eigenvalues.
"""
function EPSGetConverged(eps::SlepcEPS)
    nconv = Ref{PetscInt}(0)

    error = ccall((:EPSGetConverged, libslepc), PetscErrorCode, (CEPS, Ref{PetscInt}), eps, nconv)
    @assert iszero(error)

    return nconv[]
end

"""
    EPSGetEigenvalue(eps::SlepcEPS, ieig)

Wrapper for EPSGetEigenvalue . SLEPc 0-based indexing is used : `0 < ieig < EPSGetConverged-1` ()

A tuple (real, imag) is returned.
"""
function EPSGetEigenvalue(eps::SlepcEPS, ieig)
    eigr = Ref{PetscScalar}()
    eigi = Ref{PetscScalar}()

    error = ccall((:EPSGetEigenvalue, libslepc),
        PetscErrorCode,
        (CEPS, PetscInt, Ref{PetscScalar}, Ref{PetscScalar}),
        eps, PetscInt(ieig), eigr, eigi
    )
    @assert iszero(error)

    return eigr[], eigi[]
end

"""
EPSGetEigenvector(eps::SlepcEPS, ivec, vecr::PetscVec, veci::PetscVec)

Wrapper for `EPSGetEigenvector`. SLEPc 0-based indexing is used : `0 < ivec < EPSGetConverged-1` ()
"""
function EPSGetEigenvector(eps::SlepcEPS, ivec, vecr::PetscVec, veci::PetscVec)
    error = ccall((:EPSGetEigenvector, libslepc),
        PetscErrorCode,
        (CEPS, PetscInt, CVec, CVec),
        eps, ivec, vecr, veci
    )
    @assert iszero(error)
end

"""
    EPSGetTolerances(eps::SlepcEPS)

Wrapper for `EPSGetTolerances`
"""
function EPSGetTolerances(eps::SlepcEPS)
    tol = Ref{PetscReal}(0)
    maxits = Ref{PetscInt}(0)

    error = ccall((:EPSGetTolerances, libslepc), PetscErrorCode, (CEPS, Ref{PetscReal}, Ref{PetscInt}), eps, tol, maxits)
    @assert iszero(error)

    return tol[], maxits[]
end

"""
    EPSGetEigenpair(eps::SlepcEPS, ieig, vecr::PetscVec, veci::PetscVec)

Wrapper for `EPSGetEigenpair`. SLEPc 0-based indexing is used : `0 < ieig < EPSGetConverged-1`
"""
function EPSGetEigenpair(eps::SlepcEPS, ieig, vecr::PetscVec, veci::PetscVec)
    eigr = Ref{PetscScalar}(0.0)
    eigi = Ref{PetscScalar}(0.0)

    error = ccall((:EPSGetEigenpair, libslepc),
        PetscErrorCode,
        (CEPS, PetscInt, Ref{PetscScalar}, Ref{PetscScalar}, CVec, CVec),
        eps, PetscInt(ieig), eigr, eigi, vecr, veci
    )
    @assert iszero(error)

    return eigr[], eigi[], vecr, veci
end

"""
    EPSGetEigenpair(eps::SlepcEPS, mat::PetscMat, ieig)

Wrapper for EPSGetEigenpair without providing pre-allocated vec. but providing one of the operator matrix
"""
function EPSGetEigenpair(eps::SlepcEPS, mat::PetscMat, ieig)
    eigr = Ref{PetscScalar}()
    eigi = Ref{PetscScalar}()
    vecr, veci = MatCreateVecs(mat)

    return EPSGetEigenpair(eps, ieig, vecr, veci)
end
"""
    EPSGetEigenpair(eps::SlepcEPS, ieig)

Wrapper for EPSGetEigenpair without providing pre-allocated vec. nor one of the operator matrix
"""
function EPSGetEigenpair(eps::SlepcEPS, ieig)
    A, B = EPSGetOperators(eps)
    vecr, veci = MatCreateVecs(A)

    return EPSGetEigenpair(eps, ieig, vecr, veci)
end

"""
    EPSView(eps::SlepcEPS, viewer::PetscViewer = PetscViewerStdWorld())

Wrapper for EPSView
"""
function EPSView(eps::SlepcEPS, viewer::PetscViewer=PetscViewerStdWorld())
    error = ccall((:EPSView, libslepc), PetscErrorCode, (CEPS, CViewer), eps, viewer)
    @assert iszero(error)
end

"""
    EPSDestroy(eps::SlepcEPS)

Wrapper for EPSDestroy
"""
function EPSDestroy(eps::SlepcEPS)
    error = ccall((:EPSDestroy, libslepc), PetscErrorCode, (Ptr{CEPS},), eps.ptr)
    @assert iszero(error)
end