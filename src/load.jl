# Absolute path to libslepc.so
const libslepc = string(ENV["SLEPC_DIR"], "/", ENV["PETSC_ARCH"], "/lib/libslepc.so")