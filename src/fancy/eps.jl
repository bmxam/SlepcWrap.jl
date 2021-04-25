"""
    create_eps(A::PetscMat)

For a standard eigenvalue prolem.
"""
function create_eps(A::PetscMat; auto_setup = false)
    eps = EPSCreate()
    EPSSetOperators(eps, A)

    if (auto_setup)
        set_from_options!(eps)
        set_up!(eps)
    end

    return eps
end

"""
    create_eps(A::PetscMat, B::PetscMat)

For a generalized eigenvalue problem
"""
function create_eps(A::PetscMat, B::PetscMat; auto_setup = false)
    eps = EPSCreate()
    EPSSetOperators(eps, A, B)

    if (auto_setup)
        set_from_options!(eps)
        set_up!(eps)
    end

    return eps
end

PetscWrap.set_up!(eps::SlepcEPS) = EPSSetUp(eps)

PetscWrap.set_from_options!(eps::SlepcEPS) = EPSSetFromOptions(eps)

PetscWrap.solve!(eps::SlepcEPS) = EPSSolve(eps)

const neigs = EPSGetConverged

"""
    get_eig(eps::SlepcEPS, ieig::Integer)

Get the `ieig`-th eigenvalue as a complex number. Starts at 1 (Julia 1-based indexing)
"""
function get_eig(eps::SlepcEPS, ieig::Integer)
    eigr, eigi = EPSGetEigenvalue(eps, ieig - 1)
    return isa(eigr, Real) ? complex(eigr, eigi) : eigr
end

"""
    get_eig(eps::SlepcEPS, ieigs)

Get the `ieig`-th eigenvalues as a complex array. Starts at 1 (Julia 1-based indexing)
"""
function get_eig(eps::SlepcEPS, ieigs)
    eigs = zeros(Complex, size(ieigs))
    for i in ieigs
        eigs[i] = get_eig(eps, ieigs[i])
    end
    return eigs
end


const get_eigenvalue = get_eig

"""
    get_eigs(eps::SlepcEPS)

Returns all the converged eigenvalues as one complex array
"""
get_eigs(eps::SlepcEPS) = get_eig(eps, 1:neigs(eps))

const get_eigenvalues = get_eigs

const get_tolerances = EPSGetTolerances

get_eigenpair(eps::SlepcEPS, ieig) = EPSGetEigenpair(eps, ieig - 1)

PetscWrap.destroy!(eps::SlepcEPS) = EPSDestroy(eps)