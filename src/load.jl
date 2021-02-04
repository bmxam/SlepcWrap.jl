"""
    Function SLEPc lib location.
"""
function get_slepc_location()
    SLEPC_DIR   = haskey(ENV,"SLEPC_DIR") ? ENV["SLEPC_DIR"] : "/usr/lib/slepc"
    PETSC_ARCH  = haskey(ENV,"PETSC_ARCH") ? ENV["PETSC_ARCH"] : ""
    SLEPC_LIB = ""

    # Check SLEPC_DIR exists
    if isdir(SLEPC_DIR)

        # Define default paths
        SLEPC_LIB_DIR = joinpath(SLEPC_DIR,PETSC_ARCH,"lib")

        # Check SLEPC_LIB (.../libslepc.so or .../libslepc_real.so file) exists
        if isfile(joinpath(SLEPC_LIB_DIR,"libslepc.so"))
            SLEPC_LIB = joinpath(SLEPC_LIB_DIR,"libslepc.so")
        elseif isfile(joinpath(SLEPC_LIB_DIR,"libslepc_real.so"))
            SLEPC_LIB = joinpath(SLEPC_LIB_DIR,"libslepc_real.so")
        end
    end

    # PETSc lib not found
    if(length(SLEPC_LIB) == 0)
        # Workaround for automerging on RegistryCI
        if(haskey(ENV,"JULIA_REGISTRYCI_AUTOMERGE"))
            SLEPC_LIB = "JULIA_REGISTRYCI_AUTOMERGE"
        else
            throw(ErrorException("PETSc shared library (libslepc.so) not found. Please check that SLEPC_DIR and PETSC_ARCH env. variables are set."))
        end
    end

    return SLEPC_LIB
end


# Absolute path to libslepc.so
const libslepc = get_slepc_location()