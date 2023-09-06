"""
    create_eps(A::Mat)

For a standard eigenvalue prolem.
"""
function create_eps(A::Mat; auto_setup=false, comm::MPI.Comm=MPI.COMM_WORLD)
    eps = EPSCreate(comm)
    EPSSetOperators(eps, A)

    if (auto_setup)
        set_from_options!(eps)
        set_up!(eps)
    end

    return eps
end

"""
    create_eps(A::Mat, B::Mat)

For a generalized eigenvalue problem
"""
function create_eps(A::Mat, B::Mat; auto_setup=false, comm::MPI.Comm=MPI.COMM_WORLD)
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

const set_target! = EPSSetTarget

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
    A, _ = EPSGetOperators(eps)
    vecr, veci = MatCreateVecs(A)
    EPSGetEigenvector(eps, ivec - 1, vecr, veci)
    return vecr, veci
end

PetscWrap.destroy!(eps::SlepcEPS) = EPSDestroy(eps)

"""
    eigenvalues2file(eps::SlepcEPS, ieigs, eigs_path = "eigenvalues.dat"; mpi_rank = 0, delim=",")

Write eigenvalues to a CSV file.

If `two_cols == true`, the output file is a `(neig, 2)` array (without any complex number).

`mpi_rank` is the rank of the processor writing the eigenvalue file.
"""
function eigenvalues2file(
    eps::SlepcEPS,
    eigs_path::String="eigenvalues.dat",
    ieigs=1:neigs(eps);
    two_cols=false,
    write_index=false,
    write_header=false,
    comment="#",
    mpi_rank=0,
    delim=",")

    # Get eigenvalues
    λ = get_eig(eps, ieigs)

    # Split real and imag parts (optionnal)
    two_cols && (λ = hcat(real.(λ), imag.(λ)))

    # Add index (optionnal)
    write_index && (λ = hcat(ieigs, λ))

    # Prepare header
    if (write_header)
        header = comment
        write_index && (header *= "i,")
        if (!two_cols)
            header *= "ω"
        else
            header *= "ωr,ωi"
        end
    end

    # Write eigenvalues to file
    if (MPI.Comm_rank(eps.comm) == mpi_rank)
        open(eigs_path, "w") do io
            write_header && println(io, header)
            writedlm(io, λ, delim)
        end
    end
end

"""
    eigenvectors2file(eps::SlepcEPS, ivecs, vectors_path = "eigenvectors"; mpi_rank = 0, type = "ascii", format = PETSC_VIEWER_ASCII_CSV)

Concatenate specified eigenvectors in two files : real and imag parts.

# Warning
This function is experimental : it may allocate a lot of memory. Use it at your own risk.
"""
function eigenvectors2file(eps::SlepcEPS, vectors_path::String="eigenvectors", ivecs=1:neigs(eps); type="ascii", format=PETSC_VIEWER_ASCII_CSV)
    mat_r, mat_i = eigvecs_to_arrays(eps, ivecs)

    # Write matrices to file
    mat2file(mat_r, vectors_path * "_r.dat")
    mat2file(mat_i, vectors_path * "_i.dat")

    # Free memory
    destroy!.((mat_r, mat_i))
end

"""
    eigvecs_to_matrix(::Type{Mat}, eps::SlepcEPS, ivecs=1:neigs(eps), nrows_l = 0)
    eigvecs_to_matrix(::Type{Matrix{T}}, eps::SlepcEPS, ivecs=1:neigs(eps), nrows_l = 0) where T

Concatenate specified eigenvectors in two matrices : real and imag parts.

# Implementation
For some unknown reason, the code is crashing when trying to allocate vectors from A, B when
type of A, B is "Shell". To avoid this crash, the user can provide the number of local rows.
"""
function eigvecs_to_arrays(::Type{Mat}, eps::SlepcEPS, ivecs=1:neigs(eps), nrows_l=0)
    # Get local size
    A, B = EPSGetOperators(eps)
    irows = get_urange(A)
    nrows_l = length(irows)
    @assert nrows_l == length(get_urange(B)) # necessary condition

    # Create dense matrices (nnodes x neigs)
    mat_r = MatCreateDense(eps.comm, nrows_l, PETSC_DECIDE, PETSC_DECIDE, length(ivecs))
    mat_i = MatCreateDense(eps.comm, nrows_l, PETSC_DECIDE, PETSC_DECIDE, length(ivecs))
    set_up!.((mat_r, mat_i))

    # Allocate vectors
    vecr = create_vector(comm=eps.comm; nrows_loc=_nrows_l, nrows_glo=PETSC_DECIDE, autosetup=true)
    veci = create_vector(comm=eps.comm; nrows_loc=_nrows_l, nrows_glo=PETSC_DECIDE, autosetup=true)

    # Fill these matrices
    for (icol, ivec) in enumerate(ivecs)
        # Retrieve eigenvector (real and imag)
        #vecr, veci = get_eigenvector(eps, ieig)
        EPSGetEigenvector(eps, ivec - 1, vecr, veci)

        # Real part
        #- Convert to julia array
        array, array_ref = VecGetArray(vecr)

        #- Append to matrix -> problem using MatSetValues...
        for (iloc, iglob) in enumerate(irows)
            mat_r[iglob, icol] = array[iloc]
        end
        #MatSetValues(mat_r, irows, ivec, array, INSERT_VALUES)

        # Free memory
        VecRestoreArray(vecr, array_ref)

        # Imag part
        #- Convert to julia array
        array, array_ref = VecGetArray(veci)

        #- Append to matrix
        for (iloc, iglob) in enumerate(irows)
            mat_i[iglob, icol] = array[iloc]
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

    return mat_r, mat_i
end

function eigvecs_to_arrays(::Type{Matrix{T}}, eps::SlepcEPS, ivecs=1:neigs(eps), nrows_l=0) where {T}
    # Get local size
    if nrows_l > 0
        _nrows_l = PetscInt(nrows_l)
    else
        A, _ = EPSGetOperators(eps)
        nrows_l, ncols_l = getLocalSize(A)
        @assert nrows_l == ncols_l "Non square matrix not supported"
    end

    # Create dense matrices (ndofs x neigs)
    mat_r = zeros(T, _nrows_l, length(ivecs))
    mat_i = similar(mat_r)

    # Allocate vectors
    vecr = create_vector(eps.comm; nrows_loc=_nrows_l, nrows_glo=PETSC_DECIDE, autosetup=true)
    veci = create_vector(eps.comm; nrows_loc=_nrows_l, nrows_glo=PETSC_DECIDE, autosetup=true)

    # Fill these matrices
    for (icol, ivec) in enumerate(ivecs)
        # Retrieve eigenvector (real and imag)
        EPSGetEigenvector(eps, ivec - 1, vecr, veci)

        # Real part
        array, array_ref = getArray(vecr)
        mat_r[:, icol] .= array
        restoreArray(vecr, array_ref)

        # Imag part
        array, array_ref = getArray(veci)
        mat_i[:, icol] .= array
        restoreArray(veci, array_ref)
    end

    # Free vectors
    destroy!.((vecr, veci))

    return mat_r, mat_i
end