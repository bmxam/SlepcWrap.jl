module Complex #hide
# # Complex numbers
# Trivial example to demonstrate the capability of handling complex numbers
using PetscWrap
using SlepcWrap

# Start by checking that you're using a PETSc/SLEPc build with complex number
(PetscScalar <: Real) && error("You must configure PETSc/SLEPc with complex numbers to run this example")

# Initialize SLEPc
SlepcInitialize()

# Create matrix
A = create_matrix(4, 4; auto_setup = true)

# Get rows handled by the local processor
A_rstart, A_rend = get_range(A)

# Create diag matrix with complex numbers
J = I = A_rstart:A_rend
V = im .* I
set_values!(A, I, J, V)
assemble!(A)

# Now we set up the eigenvalue solver
eps = create_eps(A; auto_setup = true)

# Then we solve
solve!(eps)

# Show converged eivenvalues
@show get_eigenvalues(eps)

# Finally, let's free the memory
destroy!(A)
destroy!(eps)

# And call finalize when you're done
SlepcFinalize()

end #hide