using PackageCompiler
using Pkg

Pkg.activate(".")
Pkg.instantiate()

# Check if the system is running linux
if !Sys.islinux()
    "The code can only be compiled on Linux for now. Exiting..."
    exit()
end

# Set base variables
name = "dokoka"
version = "1.0b"
arch = Sys.ARCH

# Get the distro name
release_file = open("/etc/os-release", "r")
_ = readuntil(release_file, "\nID=")
distro_name = readline(release_file)
close(release_file)

# Get the distro version
release_file = open("/etc/os-release", "r")
_ = readuntil(release_file, "\nVERSION_ID=")
distro_id = readline(release_file)
close(release_file)

# Create the build
build_name = "$name-$version-$distro_name-$distro_id-$arch"
create_app(".", "build/$build_name", precompile_execution_file="scripts/Precompile.jl",
            executables=["dokoka" => "main"], force=true)
