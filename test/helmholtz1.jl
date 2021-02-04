@testset "Helmholtz 1" begin

# Only on one processor...
using PetscWrap
using SlepcWrap

# Number of mesh points and mesh step
n = 21
Δx = 1. / (n - 1)

# Initialize SLEPc
SlepcInitialize("-eps_target 0 -eps_nev 5 -st_pc_factor_shift_type NONZERO -st_type sinvert")

# Create the problem matrices, set sizes and apply "command-line" options. Note that we should
# set the number of preallocated non-zeros to increase performance.
A = MatCreate()
B = MatCreate()
MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n)
MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, n, n)
MatSetFromOptions(A)
MatSetFromOptions(B)
MatSetUp(A)
MatSetUp(B)

# Get rows handled by the local processor
A_rstart, A_rend = MatGetOwnershipRange(A)
B_rstart, B_rend = MatGetOwnershipRange(B)

# Fill matrix A  with second order derivative central scheme
for i in A_rstart:A_rend
    if(i == 1)
        A[1, 1:2] = [-2., 1] / Δx^2
    elseif (i == n)
        A[n, n-1:n] = [1., -2.] / Δx^2
    else
        A[i, i-1:i+1] = [1., -2., 1.] / Δx^2
    end
end

# Fill matrix B with identity matrix
for i in B_rstart:B_rend
    B[i,i] = -1.
end

# Set boundary conditions : u(0) = 0 and u(1) = 0. Only the processor handling the corresponding rows are playing a role here.
(A_rstart == 1) && (A[1, 1:2] = [1. 0.] )
(B_rstart == 1) && (B[1,   1] = 0.      )

(A_rend == n) && (A[n, n-1:n] = [0. 1.] )
(B_rend == n) && (B[n,     n] = 0.      )

# Assemble the matrices
MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY)
MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY)

# Now we set up the eigenvalue solver
eps = EPSCreate()
EPSSetOperators(eps, A, B)
EPSSetFromOptions(eps)
EPSSetUp(eps)

# Then we solve
EPSSolve(eps)

# And finally we can inspect the solution. Let's first get the number of converged eigenvalues:
nconv = EPSGetConverged(eps)

# Then we can get/display these eigenvalues (more precisely their square root, i.e ``\simeq \omega``)
vpr, vpi = EPSGetEigenvalue(eps, 1)
@test isapprox(vpr, π; atol = 1e-2)

# We can also play with eigen vectors. First, create two Petsc vectors to allocate memory
vecr, veci = MatCreateVecs(A)

# Then loop over the eigen pairs and retrieve eigenvectors
for ieig in 1:nconv
    vpr, vpi, vecpr, vecpi = EPSGetEigenpair(eps, ieig, vecr, veci)
end

# Finally, let's free the memory
MatDestroy(A)
MatDestroy(B)
EPSDestroy(eps)

# And call finalize when you're done
SlepcFinalize()

# Test if we reached this point
@test true

end