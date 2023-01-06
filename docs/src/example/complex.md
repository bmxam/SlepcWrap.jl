
# Complex numbers
Trivial example to demonstrate the capability of handling complex numbers

````julia
using PetscWrap
using SlepcWrap
````

Start by checking that you're using a PETSc/SLEPc build with complex number

````julia
(PetscScalar <: Real) && error("You must configure PETSc/SLEPc with complex numbers to run this example")
````

Initialize SLEPc

````julia
SlepcInitialize()
````

Create matrix

````julia
A = create_matrix(4, 4; auto_setup=true)
````

Get rows handled by the local processor

````julia
A_rstart, A_rend = get_range(A)
````

Create diag matrix with complex numbers

````julia
J = I = A_rstart:A_rend
V = im .* I
set_values!(A, I, J, V)
assemble!(A)
````

Now we set up the eigenvalue solver

````julia
eps = create_eps(A; auto_setup=true)
````

Then we solve

````julia
solve!(eps)
````

Show converged eivenvalues

````julia
@show get_eigenvalues(eps)
````

Finally, let's free the memory

````julia
destroy!(A)
destroy!(eps)
````

And call finalize when you're done

````julia
SlepcFinalize()

````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

