const CEPS = Ptr{Cvoid}

struct SlepcEPS
    ptr::Ref{CEPS}

    SlepcEPS() = new(Ref{Ptr{Cvoid}}())
end

# allows us to pass SlepcEPS objects directly into CEPS ccall signatures
Base.cconvert(::Type{CEPS}, eps::SlepcEPS) = eps.ptr[]

"""
    EPSCreate(comm, eps::SlepcEPS)

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
    create_eps(A::PetscMat)

For a standard eigenvalue prolem.
"""
function create_eps(A::PetscMat)
    eps = EPSCreate()
    EPSSetOperators(eps, A)
    return eps
end

"""
    create_eps(A::PetscMat, B::PetscMat)

For a generalized eigenvalue problem
"""
function create_eps(A::PetscMat, B::PetscMat)
    eps = EPSCreate()
    EPSSetOperators(eps, A, B)
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
PetscWrap.set_up!(eps::SlepcEPS) = EPSSetUp(eps)

"""
    EPSSetFromOptions(eps::SlepcEPS)

Wrapper for EPSSetFromOptions
"""
function EPSSetFromOptions(eps::SlepcEPS)
    error = ccall((:EPSSetFromOptions, libslepc), PetscErrorCode, (CEPS,), eps)
    @assert iszero(error)
end
PetscWrap.set_from_options!(eps::SlepcEPS) = EPSSetFromOptions(eps)

"""
    EPSGetOperators(eps::SlepcEPS)

Wrapper for EPSGetOperators
"""
function EPSGetOperators(eps::SlepcEPS)
    A = PetscMat()
    B = PetscMat()

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
PetscWrap.solve!(eps::SlepcEPS) = EPSSolve(eps)

"""
    EPSGetConverged(eps::SlepcEPS)

Wrapper for EPSGetConverged
"""
function EPSGetConverged(eps::SlepcEPS)
    nconv = Ref{PetscInt}(0)

    error = ccall((:EPSGetConverged, libslepc), PetscErrorCode, (CEPS, Ref{PetscInt}), eps, nconv)
    @assert iszero(error)

    return nconv[]
end
const neigs = EPSGetConverged

"""
    EPSGetEigenvalue(eps::SlepcEPS, ieig)

Wrapper for EPSGetEigenvalue . `ieig` must be in [1, EPSGetConverged]; it does't start at `0` as in Slepc API.

A tuple (real, imag) is returned.
"""
function EPSGetEigenvalue(eps::SlepcEPS, ieig)
    eigr = Ref{PetscScalar}()
    eigi = Ref{PetscScalar}()

    error = ccall((:EPSGetEigenvalue, libslepc),
                   PetscErrorCode,
                   (CEPS, PetscInt, Ref{PetscScalar}, Ref{PetscScalar}),
                   eps, PetscInt(ieig - 1), eigr, eigi
    )
    @assert iszero(error)

    return eigr[], eigi[]
end

"""
    get_eig(eps::SlepcEPS, ieig::Integer)

Get the `ieig`-th eigenvalue as a complex number.
"""
function get_eig(eps::SlepcEPS, ieig::Integer)
    eigr, eigi = EPSGetEigenvalue(eps, ieig)
    return complex(eigr, eigi)
end

"""
    get_eig(eps::SlepcEPS, ieigs)

Get the `ieig`-th eigenvalues as a complex array.
"""
function get_eig(eps::SlepcEPS, ieigs)
    eigs = zeros(Complex, size(ieigs))
    for i in ieigs
        eigs[i] = complex(get_eig(eps, ieigs[i]))
    end
    return eigs
end


const get_eigenvalue = get_eig

"""
    get_eigs(eps::SlepcEPS)

Returns all the converged eigenvalues as one complex array
"""
function get_eigs(eps::SlepcEPS)
    return get_eig(eps, 1:neigs(eps))
end
const get_eigenvalues = get_eigs

"""
    EPSGetTolerances(eps::SlepcEPS)

Wrapper for EPSGetTolerances
"""
function EPSGetTolerances(eps::SlepcEPS)
    tol = Ref{PetscReal}(0)
    maxits = Ref{PetscInt}(0)

    error = ccall((:EPSGetTolerances, libslepc), PetscErrorCode, (CEPS, Ref{PetscReal}, Ref{PetscInt}), eps, tol, maxits)
    @assert iszero(error)

    return tol[], maxits[]
end
const get_tolerances = EPSGetTolerances

"""
    EPSGetEigenpair(eps::SlepcEPS, ieig, vecr::PetscVec, veci::PetscVec)

Wrapper for EPSGetEigenpair. `ieig` must be in [1, EPSGetConverged]; it does't start at `0` as in Slepc API.
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
const get_eigenpair = EPSGetEigenpair

"""
    EPSView(eps::SlepcEPS, view::PetscViewer = C_NULL)

Wrapper for EPSView
"""
function EPSView(eps::SlepcEPS, viewer::PetscViewer = C_NULL)
    error = ccall((:EPSView, libslepc), PetscErrorCode, (CEPS, PetscViewer), eps, viewer)
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
PetscWrap.destroy!(eps::SlepcEPS) = EPSDestroy(eps)