
# Helmholtz equation
In this example, we use the SLEPc to find the eigenvalues of the following Helmholtz equation:
``u'' + \omega^2 u = 0`` associated to Dirichlet boundary conditions on the domain ``[0,1]``. Hence
the theoritical eigenvalues are ``\omega = k \pi`` with ``k \in \mathbb{Z}^*``; and the associated
eigenvectors are ``u(x) = \sin(k\pix)``.
A centered finite difference scheme is used for the spatial discretization.

The equation is written in matrix form ``Au = \alpha Bu`` where ``\alpha = \omega^2``.

To run this example, simplfy excute `mpirun -n your_favourite_integer julia helmholtz_FD.jl`

In this example, "fancy" names are use to interface with PETSc/SLEPc. For legacy method names check the
previous example.

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
A = create_matrix(n, n)
B = create_matrix(n, n)
set_from_options!(A)
set_from_options!(B)
set_up!(A)
set_up!(B)
```

Get rows handled by the local processor

```julia
A_rstart, A_rend = get_range(A)
B_rstart, B_rend = get_range(B)
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
assemble!(A)
assemble!(B)
```

Now we set up the eigenvalue solver

```julia
eps = create_eps(A, B)
set_from_options!(eps)
set_up!(eps)
```

Then we solve

```julia
solve!(eps)
```

And finally we can inspect the solution. Let's first get the number of converged eigenvalues:

```julia
nconv = neigs(eps)
```

Then we can get/display these eigenvalues (more precisely their square root, i.e ``\simeq \omega``)

```julia
for ieig in 1:nconv
    eig = get_eig(eps, ieig)
    @show √(real(eig))
end
```

You can also get all the converged eigenvalues in one call

```julia
eigs = get_eigenvalues(eps)
```

We can also play with eigen vectors.

```julia
for ieig in 1:nconv
    vpr, vpi, vecpr, vecpi = get_eigenpair(eps, ieig)

    # At this point, you can call VecGetArray to obtain a Julia array (see PetscWrap examples).
    # If you are on one processor, you can even plot the solution to check that you have a sinus
    # solution. On multiple processors, this would require to "gather" the solution on one processor only.
end
```

Finally, let's free the memory

```julia
destroy!(A)
destroy!(B)
destroy!(eps)
```

And call finalize when you're done

```julia
SlepcFinalize()

```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

