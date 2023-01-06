module demo #hide
# Tutorial from slepc4py manual
using PetscWrap
using SlepcWrap

# Number of points
n = 30

# Initialize SLEPc
SlepcInitialize()

# Create the problem matrices
A = MatCreate()
MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n)
MatSetFromOptions(A)
MatSetUp(A)

# Get rows handled by the local processor
rstart, rend = MatGetOwnershipRange(A)

# Fill matrix A
if rstart == 0
    MatSetValues(A, [0], [0, 1], [2.0, 1.0], INSERT_VALUES)
    rstart += 1
end
if rend == n
    MatSetValues(A, [n - 1], [n - 2, n - 1], [-1.0, 2.0], INSERT_VALUES)
    rend -= 1
end
for i in rstart:rend-1
    MatSetValues(A, [i], [i - 1, i, i + 1], [-1.0, 2.0, -1.0], INSERT_VALUES)
end

# Assemble the matrices
MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY)

# Now we set up the eigenvalue solver
eps = EPSCreate()
EPSSetOperators(eps, A)
EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE)
EPSSetFromOptions(eps)
EPSSetUp(eps)

# Then we solve
EPSSolve(eps)

# Optional : display informations aboutthe sover
EPSView(eps)

# And finally we can inspect the solution. Let's first get the number of converged eigenvalues:
nconv = EPSGetConverged(eps)

# Then we can get/display these eigenvalues
for ieig in 0:nconv-1
    vpr, vpi = EPSGetEigenvalue(eps, ieig)
    @show vpr
end

# Finally, let's free the memory
MatDestroy(A)
EPSDestroy(eps)

# And call finalize when you're done
SlepcFinalize()
end  #hide