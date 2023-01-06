```@meta
CurrentModule = SlepcWrap
```

# SlepcWrap.jl

SlepcWrap.jl is a parallel Julia wrapper for the (awesome) [SLEPc](https://slepc.upv.es/) library. As described on their main page, "SLEPc is a software library for the solution of large scale sparse eigenvalue problems on parallel computers. It is an extension of PETSc and can be used for linear eigenvalue problems in either standard or generalized form, with real or complex arithmetic.".

Note that as SLEPc is an extension of PETSc, SlepcWrap.jl is an extension of [PetscWrap.jl](https://github.com/bmxam/PetscWrap.jl).

The project is far from covering all SLEPc methods, but adding a new wrapper is very quick and easy.

Check out the examples or the API.
