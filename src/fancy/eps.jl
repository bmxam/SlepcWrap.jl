"""
    create_eps(A::PetscMat)

For a standard eigenvalue prolem.
"""
function create_eps(A::PetscMat; auto_setup = false, comm::MPI.Comm = MPI.COMM_WORLD)
    eps = EPSCreate(comm)
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
function create_eps(A::PetscMat, B::PetscMat; auto_setup = false, comm::MPI.Comm = MPI.COMM_WORLD)
    eps = EPSCreate(comm)
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

"""
    get_eigenpair(eps::SlepcEPS, ieig)

Return the `ieig`-th eigenpair (eigenvalue, eigenvector). `ieig` starts at 1 (Julia 1-based indexing)
"""
get_eigenpair(eps::SlepcEPS, ieig) = EPSGetEigenpair(eps, ieig - 1)

"""
    get_eigenvector(eps::SlepcEPS, ivec)

Return a tuple (vecr, veci) corresponding to the `ivec`-th eigenvector. `ivec` starts at 1 (Julia 1-based indexing)
"""
function get_eigenvector(eps::SlepcEPS, ivec)
    A, B = EPSGetOperators(eps)
    vecr, veci = MatCreateVecs(A)
    EPSGetEigenvector(eps, ivec - 1, vecr, veci)
    return vecr, veci
end

PetscWrap.destroy!(eps::SlepcEPS) = EPSDestroy(eps)

"""
    eigenvalues2file(eps::SlepcEPS, ieigs, eigs_path = "eigenvalues.dat"; mpi_rank = 0, delim=",")

Write eigenvalues to a CSV file.

`mpi_rank` is the rank of the processor writing the eigenvalue file.
"""
function eigenvalues2file(eps::SlepcEPS, ieigs, eigs_path = "eigenvalues.dat"; mpi_rank = 0, delim=",")
    # Write eigenvalues to file
    λ = get_eigs(eps)
    if (MPI.Comm_rank(eps.comm) == mpi_rank)
        open(eigs_path, "w") do io
            writedlm(io, λ, delim)
        end
    end
end

"""
    eigenvectors2file(eps::SlepcEPS, ivecs, vectors_path = "eigenvectors"; mpi_rank = 0, type = "ascii", format = PETSC_VIEWER_ASCII_CSV)

Concatenate specified eigenvectors in two files : real and imag parts.

`mpi_rank` is the rank of the processor writing the eigenvalue file.

# Warning
This experimental : it may allocate a lot of memory. Use it at your own risk.
"""
function eigenvectors2file(eps::SlepcEPS, ivecs, vectors_path = "eigenvectors"; mpi_rank = 0, type = "ascii", format = PETSC_VIEWER_ASCII_CSV)
    # Get local size
    A, B = EPSGetOperators(eps)
    irows = get_urange(A)
    nrows_l = length(irows)
    @assert nrows_l == length(get_urange(B)) # necessary condition

    # Create dense matrices (nnodes x neigs)
    mat_r = MatCreateDense(eps.comm, nrows_l, length(ivecs), PETSC_DECIDE, PETSC_DECIDE)
    mat_i = MatCreateDense(eps.comm, nrows_l, length(ivecs), PETSC_DECIDE, PETSC_DECIDE)
    set_up!.((mat_r, mat_i))

    # Allocate vectors
    vecr = create_vector(PETSC_DECIDE, nrows_l; comm = eps.comm)
    veci = create_vector(PETSC_DECIDE, nrows_l; comm = eps.comm)
    set_up!.((vecr, veci))

    # Fill these matrices
    for ivec in ivecs
        # Retrieve eigenvector (real and imag)
        #vecr, veci = get_eigenvector(eps, ieig)
        EPSGetEigenvector(eps, ivec - 1, vecr, veci)

        # Real part
        #- Convert to julia array
        array, array_ref = VecGetArray(vecr)

        #- Append to matrix -> problem using MatSetValues...
        for (iloc, iglob) in enumerate(irows)
            mat_r[iglob,ivec] = array[iloc]
        end
        #MatSetValues(mat_r, irows, ivec, array, INSERT_VALUES)

        # Free memory
        VecRestoreArray(vecr, array_ref)

        # Imag part
        #- Convert to julia array
        array, array_ref = VecGetArray(veci)

        #- Append to matrix
        for (iloc, iglob) in enumerate(irows)
            mat_i[iglob,ivec] = array[iloc]
        end
        #MatSetValues(mat_i, irows, ivec, array, INSERT_VALUES)

        # Free memory
        VecRestoreArray(veci, array_ref)
    end

    # Free vectors
    destroy!.((vecr, veci))

    # Assemble matrices
    assemble!(mat_r, MAT_FINAL_ASSEMBLY)
    assemble!(mat_i, MAT_FINAL_ASSEMBLY)

    # Prepare viewer
    viewer = PetscViewerCreate(eps.comm)
    set_type!(viewer, type)
    set_mode!(viewer, FILE_MODE_WRITE)
    push_format!(viewer, format)

    # Write real matrix to file
    set_name!(viewer, vectors_path * "_r.dat")
    MatView(mat_r, viewer)

    # Write imag matrix to file (reuse viewer)
    set_name!(viewer, vectors_path * "_i.dat")
    MatView(mat_i, viewer)
    destroy!(viewer)

    # Free memory
    destroy!.((mat_r, mat_i))
end