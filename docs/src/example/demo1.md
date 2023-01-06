
# Demo 1
Tutorial from slepc4py manual

````julia
using PetscWrap
using SlepcWrap
````

Number of points

````julia
n = 30
````

Initialize SLEPc

````julia
SlepcInitialize()
````

Create the problem matrices

````julia
A = MatCreate()
MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n)
MatSetFromOptions(A)
MatSetUp(A)
````

Get rows handled by the local processor

````julia
rstart, rend = MatGetOwnershipRange(A)
````

Fill matrix A

````julia
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
````

Assemble the matrices

````julia
MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY)
````

Now we set up the eigenvalue solver : note that contrary to the original
tutorial, we select the "smallest magnitude" eigenvalues

````julia
eps = EPSCreate()
EPSSetOperators(eps, A)
EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE)
EPSSetFromOptions(eps)
EPSSetUp(eps)
````

Then we solve

````julia
EPSSolve(eps)
````

Optional : display informations aboutthe sover

````julia
EPSView(eps)
````

And finally we can inspect the solution. Let's first get the number of converged eigenvalues:

````julia
nconv = EPSGetConverged(eps)
````

Then we can get/display these eigenvalues

````julia
for ieig in 0:nconv-1
    vpr, vpi = EPSGetEigenvalue(eps, ieig)
    @show vpr
end
````

Finally, let's free the memory

````julia
MatDestroy(A)
EPSDestroy(eps)
````

And call finalize when you're done

````julia
SlepcFinalize()
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

