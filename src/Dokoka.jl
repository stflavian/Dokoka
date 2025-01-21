#!/usr/bin/env -S julia --project=~/Projects/Dokoka/.

module Dokoka

export Atom, Molecule, Configuration
export parseargs, readmol
export goangular!, goradial!, gorotational!, gorandom!

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
    settings.description = "A program that generates conformations of molecules"
    settings.version = "1.0a"
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
    end
    
    if interactive 
        return ArgParse.parse_args(settings)
    end
    
    return 0
end


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


function getsize(molecule::Molecule)
    x = diff(collect(extrema([atom.r[1] for atom in molecule.atoms])))[1]
    y = diff(collect(extrema([atom.r[2] for atom in molecule.atoms])))[1]
    z = diff(collect(extrema([atom.r[3] for atom in molecule.atoms])))[1]
    return LA.norm([x, y, z])
end


function getcenter(molecule::Molecule)
    return sum([atom.r for atom in molecule.atoms]) ./ length(molecule.atoms)
end


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


function center!(molecule::Molecule)
    cm = getcenter(molecule)
    for atom in molecule.atoms
        atom.r -= cm
    end
end


function displace!(molecule::Molecule, displacement::Vector{Float64})
    for atom in molecule.atoms
        atom.r += displacement
    end
end


function goradial!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
    axis::Vector{Float64}, angle::Float64, n::Int64)
    
    # Center the molecules
    center!(molecule1)
    center!(molecule2)
    
    # Move second molecule to initial position and rotate it
    rotate!(molecule2, angle * pi / 180, axis[1], axis[2], axis[3])
    displace!(molecule2, position) 
    
    # Get unit vector and distances
    init_distance = LA.norm(position) 
    unit = position ./ init_distance
    d_distance = getsize(molecule1) / n

    # Generate the conformations and write the results
    io = open("radial.xyz", "w")
    for i in 1:n
        displace!(molecule2, unit .* d_distance)
        write_movie(Configuration([molecule1, molecule2]), io, init_distance + i * d_distance)
    end
    close(io)
end


function gorotational!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
    axis::Vector{Float64}, angle::Float64, n::Int64)
    
    # Center the molecules
    center!(molecule1)
    center!(molecule2)
    
    # Move second molecule to initial position and rotate it
    rotate!(molecule2, angle * pi / 180, axis[1], axis[2], axis[3])
    displace!(molecule2, position) 

    # Get unit of rotational movement
    d_angle = 2 * pi / n

    # Generate the conformations and write the results
    io = open("rotational.xyz", "w")
    for i in 1:n
        center!(molecule2)
        rotate!(molecule2, d_angle, axis[1], axis[2], axis[3])
        displace!(molecule2, position)
        write_movie(Configuration([molecule1, molecule2]), io, i * d_angle)
    end
    close(io)
end


function goangular!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
    axis::Vector{Float64}, angle::Float64, n::Int64)
    
    # Center the molecules
    center!(molecule1)
    center!(molecule2)
    
    # Move second molecule to initial position and rotate it
    rotate!(molecule2, angle * pi / 180, axis[1], axis[2], axis[3])
    displace!(molecule2, position) 

    # Get unit of rotational movement
    d_angle = 2 * pi / n

    # Generate the conformations and write the results
    io = open("angular.xyz", "w")
    for i in 1:n
        rotate!(molecule2, d_angle, axis[1], axis[2], axis[3])
        write_movie(Configuration([molecule1, molecule2]), io, i * d_angle)
    end
    close(io)
end


function gorandom!(molecule1::Molecule, molecule2::Molecule, position::Vector{Float64},
    axis::Vector{Float64}, angle::Float64, n::Int64)
    
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
        write_movie(Configuration([molecule1, molecule2]), io, "Generated by dokoka")
    end
    close(io)
end


function write_movie(configuration::Configuration, io::IO, comment::Any)
    size = length(configuration.molecules) * length(configuration.molecules[1].atoms)
    println(io, size)
    println(io, comment)
    for molecule in configuration.molecules
        for atom in molecule.atoms
            println(io, atom.species, " ", atom.r[1], " ", atom.r[2], " ", atom.r[3])
        end
    end
end


function julia_main()::Cint
    
    # Parse the CLI arguments
    args = parseargs()
    
    # Assign the values to the variables
    molecule1 = readmol(args["molecules"][1])
    molecule2 = readmol(args["molecules"][2])
    position = args["position"]
    axis = args["axis"]
    angle = args["angle"]
    n = args["number"]
    
    methods = Dict{String, Function}(
        "radial" => goradial!,
        "rotational" => gorotational!,
        "angular" => goangular!,
        "random" => gorandom!
    )
    
    get(methods, args["method"], gorandom!)(molecule1, molecule2, position, axis, angle, n)
    return 0
end

end


function (@main)(_)
    Dokoka.julia_main() 
end
