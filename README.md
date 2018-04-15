# Nimbus
Lightweight (restricted) Hartree-Fock program to calculate the wavefunction of the molecular orbitals. The output of the program are a series of density files which can be converted to wavefront object using the [den2obj](https://github.com/ifilot/den2obj) tool.

# Compilation instructions

Nimbus has the following dependencies:
* Boost
* Eigen3
* TCLAP

To compile, run the following commands:
```
https://github.com/ifilot/nimbus.git
cd den2obj
mkdir build
cd build
cmake ../src
make -j5
```

# Usage

Note that nimbus needs to be run from its build folder (see compilation instructions) in order to find the basis set files.

```
./nimbus -i ../molecules/co.in
```

Optional settings
```
-r   gridpoint resolution for the wavefunction amplitude files (default = 0.1)
-b   unit cell resolution (default = 10.0)
```

For more information
```
./nimbus --help

```
