[![](https://img.shields.io/badge/docs-stable-red.svg)](https://bmxam.github.io/SlepcWrap.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://bmxam.github.io/SlepcWrap.jl/dev/)

# SlepcWrap.jl

SlepcWrap.jl is a parallel Julia wrapper for the (awesome) [SLEPc](https://slepc.upv.es/) library. As described on their main page, "SLEPc is a software library for the solution of large scale sparse eigenvalue problems on parallel computers. It is an extension of PETSc and can be used for linear eigenvalue problems in either standard or generalized form, with real or complex arithmetic.".

Note that as SLEPc is an extension of PETSc, SlepcWrap.jl is an extension of [PetscWrap.jl](https://github.com/bmxam/PetscWrap.jl).

The project is far from covering all SLEPc methods, but adding a new wrapper is very quick and easy.
## How to install it
You must have installed the SLEPc library (and necessarily the PETSc library as wall) on your computer and set the two following environment variables : `SLEPC_DIR` and `PETSC_ARCH`.

At run time, PetscWrap.jl looks for the `libslepc.so` using these environment variables and "load" the library.

To install the package, use the Julia package manager:
```Julia
pkg> add SlepcWrap
```
## Contribute
Any contribution(s) and/or remark(s) are welcome! If you need a function that is not wrapped yet but you don't think you are capable of contributing, post an issue with a minimum working example.

## SLEPc compat.
This version of PetscWrap.jl has been tested with slepc-3.13.1 Complex numbers are not supported yet.

## How to use it
You will find examples of use by building the documentation: `julia SlepcWrap.jl/docs/make.jl`. Here is one of the examples:
### Helmholtz equation
In this example, we use the SLEPc to find the eigenvalues of the following Helmholtz equation:
``u'' + \omega^2 u = 0`` associated to Dirichlet boundary conditions on the domain ``[0,1]``. Hence
the theoritical eigenvalues are ``\omega = k \pi`` with ``k \in \mathbb{Z}^*``; and the associated
eigenvectors are ``u(x) = \sin(k\pix)``.
A centered finite difference scheme is used for the spatial discretization.

The equation is written in matrix form ``Au = \alpha Bu`` where ``\alpha = \omega^2``.

To run this example, simplfy excute `mpirun -n your_favourite_integer julia helmholtz_FD.jl`

In this example, PETSc/SLEPc legacy method names are used. For more fancy names, check the next example.

Note that the way we achieve things in the document can be highly improved and the purpose of this example
is only demonstrate some method calls to give an overview.

Start by importing both `PetscWrap`, for the distributed matrices, and `SlepcWrap` for the eigenvalues.

```julia
using PetscWrap
using SlepcWrap
```

Number of mesh points and mesh step

```julia
n = 21
Δx = 1. / (n - 1)
```

Initialize SLEPc. Either without arguments, calling `SlepcInitialize()` or using "command-line" arguments.
To do so, either provide the arguments as one string, for instance
`SlepcInitialize("-eps_max_it 100 -eps_tol 1e-5")` or provide each argument in
separate strings : `PetscInitialize(["-eps_max_it", "100", "-eps_tol", "1e-5")`.
Here we ask for the five closest eigenvalues to ``0``, using a non-zero pivot for the LU factorization and a
"shift-inverse" process.

```julia
SlepcInitialize("-eps_target 0 -eps_nev 5 -st_pc_factor_shift_type NONZERO -st_type sinvert")
```

Create the problem matrices, set sizes and apply "command-line" options. Note that we should
set the number of preallocated non-zeros to increase performance.

```julia
A = MatCreate()
B = MatCreate()
MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n)
MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, n, n)
MatSetFromOptions(A)
MatSetFromOptions(B)
MatSetUp(A)
MatSetUp(B)
```

Get rows handled by the local processor

```julia
A_rstart, A_rend = MatGetOwnershipRange(A)
B_rstart, B_rend = MatGetOwnershipRange(B)
```

Fill matrix A  with second order derivative central scheme

```julia
for i in A_rstart:A_rend
    if(i == 1)
        A[1, 1:2] = [-2., 1] / Δx^2
    elseif (i == n)
        A[n, n-1:n] = [1., -2.] / Δx^2
    else
        A[i, i-1:i+1] = [1., -2., 1.] / Δx^2
    end
end
```

Fill matrix B with identity matrix

```julia
for i in B_rstart:B_rend
    B[i,i] = -1.
end
```

Set boundary conditions : u(0) = 0 and u(1) = 0. Only the processor handling the corresponding rows are playing a role here.

```julia
(A_rstart == 1) && (A[1, 1:2] = [1. 0.] )
(B_rstart == 1) && (B[1,   1] = 0.      )

(A_rend == n) && (A[n, n-1:n] = [0. 1.] )
(B_rend == n) && (B[n,     n] = 0.      )
```

Assemble the matrices

```julia
MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY)
MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY)
```

Now we set up the eigenvalue solver

```julia
eps = EPSCreate()
EPSSetOperators(eps, A, B)
EPSSetFromOptions(eps)
EPSSetUp(eps)
```

Then we solve

```julia
EPSSolve(eps)
```

And finally we can inspect the solution. Let's first get the number of converged eigenvalues:

```julia
nconv = EPSGetConverged(eps)
```

Then we can get/display these eigenvalues (more precisely their square root, i.e ``\simeq \omega``)

```julia
for ieig in 1:nconv
    vpr, vpi = EPSGetEigenvalue(eps, ieig)
    @show √(vpr), √(vpi)
end
```

We can also play with eigen vectors. First, create two Petsc vectors to allocate memory

```julia
vecr, veci = MatCreateVecs(A)
```

Then loop over the eigen pairs and retrieve eigenvectors

```julia
for ieig in 1:nconv
    vpr, vpi, vecpr, vecpi = EPSGetEigenpair(eps, ieig, vecr, veci)

    # At this point, you can call VecGetArray to obtain a Julia array (see PetscWrap examples).
    # If you are on one processor, you can even plot the solution to check that you have a sinus
    # solution. On multiple processors, this would require to "gather" the solution on one processor only.
end
```

Finally, let's free the memory

```julia
MatDestroy(A)
MatDestroy(B)
EPSDestroy(eps)
```

And call finalize when you're done

```julia
SlepcFinalize()

```

