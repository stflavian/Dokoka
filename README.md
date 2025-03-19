# Dokoka (何処か)

A Julia script for generating conformations of ionic liquids using random positioning 
of cations and anions.

## Table of Contents

1. [Introduction](#introduction)
2. [Usage](#usage)
3. [Installation](#installation)
4. [Contributing](#contributing)

## Introduction
-------------

Dokoka is a script that generates conformations of ionic liquids using random or 
user-defined positioning of cations and anions. The script uses Julia, a high-level, 
high-performance language, and is compiled to ensure maximal performance. 

## Usage
-----

To use the script, simply run it from the command line with the following arguments:

*   `molecules` : The path to the structure of the input molecules, usually the 
cation followed by the anion (eg. `cation.xyz ../folder/anion.xyz`)
*   `-m` or `--method` : The method of generation; radial, rotational, angular or 
random
*   `-r` or `--position` : The starting position of the anion with respect to the 
cation in vector format (eg. `1.4 2.5 -0.4`)
*   `-x` or `--axis` : The axis of rotation (eg. `0 0 1`)
*   `-a` or `--angle` : The angle of rotation in degrees (eg. `180`)
*   `-l` or `--limit` : The movement limit (distance for radial motion, angle for
angular and rotational motion) 
*   `-n` or `--number`: The number of conformations generated

Example:

```bash
dokoka C6H12O6.xyz C7H5NO2.xyz -m radial -r 2.3 3.1 -0.7 --number 1000
```

```bash
dokoka C4MIM.xyz Cl.xyz --method random --position 0.1 0.0 0.0 --number 1000
```

## Installation
-------------

Pre-compiled binaries are included in the [Release section][1] on GitHub. You can 
download and use these binaries as is. For maximal compatibility, use the Debian 12
binary.

To build the script manually, follow these steps:

1. Clone this repository using `git`.
2. Navigate to the project directory.
3. Run `julia scripts/CreateImage.jl` to compile the binaries.

The compiled scripts will be available in the `build/{build-name}/bin/` directory.

## Contributing
------------

Feel free to contribute to this script by submitting a pull request or creating an 
issue.

[1]: https://github.com/stflavian/Dokoka/releases

