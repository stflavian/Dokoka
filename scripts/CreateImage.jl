using PackageCompiler
using Pkg

Pkg.activate(".")
Pkg.instantiate()

if !Sys.islinux()
    "The code can only be compiled on Linux for now. Exiting..."
    exit()
end

release_file = open("/etc/os-release", "r")
_ = readuntil(release_file, "\nID=")
distro_name = readline(release_file)
close(release_file)

create_app(".", "build/$distro_name", precompile_execution_file="scripts/Precompile.jl",
            executables=["dokoka" => "julia_main"])
