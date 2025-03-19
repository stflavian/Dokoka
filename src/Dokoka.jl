#!/usr/bin/env -S julia --project=~/Projects/Dokoka/.

module Dokoka

export parseargs
export goangular!, goradial!, gorotational!, gorandom!
export julia_main

import LinearAlgebra as LA
import ArgParse
import Random

mutable struct Atom
    species::String
    r::Vector{Float64}
end

mutable struct Molecule
    atoms::Vector{Atom}
end

mutable struct Configuration
    molecules::Vector{Molecule}
end


function parseargs(interactive::Bool = true)
    settings = ArgParse.ArgParseSettings()
    settings.prog = "dokoka"
    settings.description = "A program that generates conformations of ionic liquids"
    settings.version = "1.0b"
    settings.add_version = true

    ArgParse.@add_arg_table! settings begin
        "molecules"
            required = true
            nargs = 2
            action = "store_arg"
            arg_type = String
            help = "Molecules to be repositioned"
        "--method", "-m"
            action = "store_arg"
            arg_type = String
            default = "radial"
            help = "Type of displacement used for the second molecule"
        "--position", "-r"
            nargs = 3
            action = "store_arg"
            arg_type = Float64
            default = [0.0, 0.0, 0.0]
            help = "Position of the second molecule"
        "--axis", "-x"
            nargs = 3
            action = "store_arg"
            arg_type = Float64
            default = [0.0, 0.0, 1.0]
            help = "Axis of rotation used for the second molecule"
        "--angle", "-a"
            action = "store_arg"
            arg_type = Float64
            default = 0.0
            help = "Angle of rotation used for the second molecule"
        "--number", "-n"
            action = "store_arg"
            arg_type = Int64
            default = 10
            help = "Number of conformations generated"
        "--limit", "-l"
            action = "store_arg"
            arg_type = Float64
            default = 0
            help = "The upper limit for the movement (ignored for random motion)"
    end
    
    if interactive 
        return ArgParse.parse_args(settings)
    end
    
    return 0
end


"""
    readmol(filename::String)

Read an .xyz file and output a Molecule object from it.
"""
function readmol(filename::String)
    io = open(filename, "r")
    natoms = parse(Int64, readline(io))
    atoms = Vector{Atom}(undef, natoms)
    readline(io)
    for i in 1:natoms
        data = split(readline(io))
        atoms[i] = Atom(data[1], [parse(Float64, data[2]); parse(Float64, data[3]); 
                                      parse(Float64, data[4])])
    end
    close(io)
    return Molecule(atoms)
end


"""
    getsize(molecule::Molecule)

Enclose the molecule in a box and get the long diagonal.
"""
function getsize(molecule::Molecule)
    x = diff(collect(extrema([atom.r[1] for atom in molecule.atoms])))[1]
    y = diff(collect(extrema([atom.r[2] for atom in molecule.atoms])))[1]
    z = diff(collect(extrema([atom.r[3] for atom in molecule.atoms])))[1]
    return LA.norm([x, y, z])
end


"""
    getcenter(molecule::Molecule)

Get the center of the molecule using unitary mass.
"""
function getcenter(molecule::Molecule)
    return sum([atom.r for atom in molecule.atoms]) ./ length(molecule.atoms)
end


"""
    rotate!(molecule::Molecule, angle::Float64, i::Float64, j::Float64, k::Float64)

Rotate the molecule around an axis which passes through the origin.
"""
function rotate!(molecule::Molecule, angle::Float64, i::Float64, j::Float64, k::Float64)
    
    rotation_vector = sin(angle/2) / sqrt(i^2 + j^2 + k^2) * [i; j; k]
    q0 = cos(angle/2)
    qlen2 = rotation_vector[1]^2 + rotation_vector[2]^2 + rotation_vector[3]^2
    
    for a = eachindex(molecule.atoms)
        molecule.atoms[a].r = (q0^2 - qlen2) * molecule.atoms[a].r + 
                    2 * LA.dot(rotation_vector, molecule.atoms[a].r) * rotation_vector + 
                    2 * q0 * LA.cross(rotation_vector, molecule.atoms[a].r)
    end
end


"""
    center!(molecule::Molecule)

Place the center of the molecule in the origin.

See also [`getcenter`](@ref).
"""
function center!(molecule::Molecule)
    cm = getcenter(molecule)
    for atom in molecule.atoms
        atom.r -= cm
    end
end


"""
    displace!(molecule::Molecule)

Displace a molecule by a given vector.
"""
function displace!(molecule::Molecule, displacement::Vector{Float64})
    for atom in molecule.atoms
        atom.r += displacement
    end
end


"""
    check_overlap(molecule1::Molecule, molecule2::Molecule)

Check if the atoms of the two molecules overlap.
"""
function check_overlap(molecule1::Molecule, molecule2::Molecule)
    for atom1 in molecule1.atoms
        for atom2 in molecule2.atoms
            if LA.norm(atom1.r - atom2.r) <= 2.5
                return true
            end
        end
    end
    return false
end


"""
    untangle!(molecule1::Molecule, molecule2::Molecule)

Untangle the two molecules if they overlap.

Move the second molecule by 0.1 Angstroms along the axis connecting the centers of the two 
molecules until the two molecules no longer overlap.
"""
function untangle!(molecule1::Molecule, molecule2::Molecule)
    center1 = getcenter(molecule1)
    center2 = getcenter(molecule2)
    
    distance = LA.norm(center1 - center2)
    unit = (center2 - center1) ./ distance
    
    while check_overlap(molecule1, molecule2)
        displace!(molecule2, unit * 0.1)
    end
end


"""
    goradial!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
              axis::Vector{Float64}, angle::Float64, n::Int64, limit::Float64)

Radially displace the second molecule along the axis connecting the centers of the two 
molecules. 

The second molecule can be rotated alongside its center using the `axis` and `angle` 
keywords. Since the user has full control over the initial position, untangling is not
used. The molecule will be displaced until the distance between the two centers is larger
than 10 Angstroms, or the user defined limit.
"""
function goradial!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
    axis::Vector{Float64}, angle::Float64, n::Int64, limit::Float64)
    
    # Set limit
    if limit == 0
        limit = 10
    end

    # Center the molecules
    center!(molecule1)
    center!(molecule2)
    
    # Move second molecule to initial position and rotate it
    rotate!(molecule2, angle * pi / 180, axis[1], axis[2], axis[3])
    displace!(molecule2, position) 
    
    # Get unit vector and distances
    init_distance = LA.norm(position) 
    unit = position ./ init_distance
    d_distance = (limit - init_distance) / n

    # Generate the conformations and write the results
    io = open("radial.xyz", "w")
    for i in 0:n-1
        displace!(molecule2, unit .* d_distance)
        write_movie(Configuration([molecule1, molecule2]), io, init_distance + i * d_distance)
    end
    close(io)
end


"""
    gorotational!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
                  axis::Vector{Float64}, angle::Float64, n::Int64, limit::Float64)

Rotate the second molecule in place along a user defined axis. 

The second molecule can be rotated alongside its center using the `axis` and `angle` 
keywords. Since the user has full control over the initial position, untangling is not
used. The molecule will be rotated 360 degrees, or the user defined limit.
"""
function gorotational!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
    axis::Vector{Float64}, angle::Float64, n::Int64, limit::Float64)
    
    # Set limit
    if limit == 0
        limit = 360
    end

    # Center the molecules
    center!(molecule1)
    center!(molecule2)
    
    # Move second molecule to initial position and rotate it
    rotate!(molecule2, angle * pi / 180, axis[1], axis[2], axis[3])
    displace!(molecule2, position) 

    # Get unit of rotational movement
    d_angle = limit * pi / 180 / n

    # Generate the conformations and write the results
    io = open("rotational.xyz", "w")
    for i in 0:n-1
        center!(molecule2)
        rotate!(molecule2, d_angle, axis[1], axis[2], axis[3])
        displace!(molecule2, position)
        write_movie(Configuration([molecule1, molecule2]), io, i * d_angle)
    end
    close(io)
end


"""
    goangular!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
               axis::Vector{Float64}, angle::Float64, n::Int64, limit::Float64)

Rotate the second molecule around the origin along a user defined axis. 

The second molecule can be rotated alongside its center using the `axis` and `angle` 
keywords. Untangling is used. The molecule will be rotated 360 degrees, or the user defined 
limit.
"""
function goangular!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
    axis::Vector{Float64}, angle::Float64, n::Int64, limit::Float64)
    
    # Set limit
    if limit == 0
        limit = 360
    end

    # Center the molecules
    center!(molecule1)
    center!(molecule2)
    
    # Move second molecule to initial position and rotate it
    rotate!(molecule2, angle * pi / 180, axis[1], axis[2], axis[3])
    displace!(molecule2, position) 

    # Get unit of rotational movement
    d_angle = limit * pi / 180 / n

    # Generate the conformations and write the results
    io = open("angular.xyz", "w")
    for i in 0:n-1
        rotate!(molecule2, d_angle, axis[1], axis[2], axis[3])
        untangle!(molecule1, molecule2)
        write_movie(Configuration([molecule1, molecule2]), io, i * d_angle)
    end
    close(io)
end


"""
    gorandom!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
              axis::Vector{Float64}, angle::Float64, n::Int64, limit::Float64)

Rotate the second molecule randomly around the first molecule. 
"""
function gorandom!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
    axis::Vector{Float64}, angle::Float64, n::Int64, limit::Float64)

    # Center the molecules
    center!(molecule1)
    center!(molecule2)
    
    # Move second molecule to initial position and rotate it
    rotate!(molecule2, angle * pi / 180, axis[1], axis[2], axis[3])
    displace!(molecule2, position) 
    
    is = 2 * rand(Float64, n) .- 1
    js = 2 * rand(Float64, n) .- 1
    ks = 2 * rand(Float64, n) .- 1
    angles = 2 * pi * rand(Float64, n)

    io = open("random.xyz", "w")
    for i = 1:n
        rotate!(molecule2, angles[i], is[i], js[i], ks[i])
        untangle!(molecule1, molecule2)
        write_movie(Configuration([molecule1, molecule2]), io, "Generated by dokoka")
        center!(molecule2)
        displace!(molecule2, position)
    end
    close(io)
end


"""
    write_movie(configuration::Configuration, io::IO, comment::Any)

Append the conformation to the specified IO file.
"""
function write_movie(configuration::Configuration, io::IO, comment::Any)
    size = length(configuration.molecules[1].atoms) + length(configuration.molecules[2].atoms)
    println(io, size)
    println(io, comment)
    for molecule in configuration.molecules
        for atom in molecule.atoms
            println(io, atom.species, " ", atom.r[1], " ", atom.r[2], " ", atom.r[3])
        end
    end
end

"""
    main()::Cint

Julia function used for compilation.
"""
function main()::Cint
    
    # Parse the CLI arguments
    args = parseargs()
    
    # Assign the values to the variables
    molecule1 = readmol(args["molecules"][1])
    molecule2 = readmol(args["molecules"][2])
    position = args["position"]
    axis = args["axis"]
    angle = args["angle"]
    n = args["number"]
    limit = args["limit"]
    
    methods = Dict{String, Function}(
        "radial" => goradial!,
        "rotational" => gorotational!,
        "angular" => goangular!,
        "random" => gorandom!
    )
    
    get(methods, args["method"], gorandom!)(molecule1, molecule2, position, 
                                            axis, angle, n, limit)
    return 0
end

#main()

end
