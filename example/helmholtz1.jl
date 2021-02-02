module Helmholtz1 #hide
# # Helmholtz equation
# In this example, we use the SLEPc to find the eigenvalues of the following Helmholtz equation:
# ``u'' + \omega^2 u = 0`` associated to Dirichlet boundary conditions on the domain ``[0,1]``. Hence
# the theoritical eigenvalues are ``\omega = k \pi`` with ``k \in \mathbb{Z}^*``; and the associated
# eigenvectors are ``u(x) = \sin(k\pix)``.
# A centered finite difference scheme is used for the spatial discretization.
#
# The equation is written in matrix form ``Au = \alpha Bu`` where ``\alpha = \omega^2``.
#
# To run this example, simplfy excute `mpirun -n your_favourite_integer julia helmholtz_FD.jl`
#
# In this example, PETSc/SLEPc legacy method names are used. For more fancy names, check the next example.
#
# Note that the way we achieve things in the document can be highly improved and the purpose of this example
# is only demonstrate some method calls to give an overview.
#
# Start by importing both `PetscWrap`, for the distributed matrices, and `SlepcWrap`
using PetscWrap
using SlepcWrap

# Number of mesh points and mesh step
n = 21
Δx = 1. / (n - 1)

# Initialize SLEPc. Either without arguments, calling `SlepcInitialize()` or using "command-line" arguments.
# To do so, either provide the arguments as one string, for instance
# `SlepcInitialize("-eps_max_it 100 -eps_tol 1e-5")` or provide each argument in
# separate strings : `PetscInitialize(["-eps_max_it", "100", "-eps_tol", "1e-5")`.
# Here we ask for the five closest eigenvalues to ``0``, using a non-zero pivot for the LU facorization.
SlepcInitialize("-eps_target 0 -eps_nev 5 -st_pc_factor_shift_type NONZERO ")

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
@show nconv

# Then we can get/display these eigenvalues (more precisely their square root, i.e ``\simeq \omega``)
for ieig in 1:nconv
    vpr, vpi = EPSGetEigenvalue(eps, ieig)
    @show √(vpr), √(vpi)
end

# vecr, veci = MatCreateVecs(A)

# for ieig in 1:nconv
#     println("Eigenvalue $ieig")
#     vpr, vpi = EPSGetEigenvalue(eps, ieig)
#     @show vpr, vpi

#     vpr, vpi, vecpr, vecpi = EPSGetEigenpair(eps, ieig, vecr, veci)
#     @show vpr, vpi

#     vpr, vpi, vecpr, vecpi = EPSGetEigenpair(eps, A, ieig)
#     @show vpr, vpi
# end

# Finally, let's free the memory
MatDestroy(A)
MatDestroy(B)
EPSDestroy(eps)

# And call finalize when you're done
SlepcFinalize()
end #hide