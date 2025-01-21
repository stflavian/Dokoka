using PackageCompiler
using Pkg

Pkg.activate(".")
Pkg.instantiate()

create_app("Dokoka", "../build", filter_stdlibs=true, precompile_execution_file="scripts/precompile.jl")
