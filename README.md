# Nimbus
Lightweight (restricted) Hartree-Fock program to calculate the wavefunction of the molecular orbitals. The output of the program are a series of density files which can be converted to wavefront object using the [den2obj](https://github.com/ifilot/den2bin) tool.

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

Example output:
```
--------------------------------------------------------------
Executing nimbus v.0.2.0
Author: Ivo Filot <i.a.w.filot@tue.nl>
--------------------------------------------------------------

           Reading input file
========================================
Reading file:       ../molecules/co.in

Atoms in system:    2
========================================
6      0.000000      0.000000     -1.058000
8      0.000000      0.000000      1.058000
========================================
Total number of GTOs: 44

-------------------------------------------------
   1 |  22.68431001
   2 | -62.05768928
   3 | -89.93896303
   4 | -101.60365294
   5 | -107.14224778
   6 | -109.89170878
   7 | -111.27285405
   8 | -111.96744659
   9 | -112.31632510
  10 | -112.49132327
  11 | -112.57901573
  12 | -112.62292991
  13 | -112.64491168
  14 | -112.65591187
  15 | -112.66141562
  16 | -112.66416898
  17 | -112.66554629
  18 | -112.66623521
  19 | -112.66657980
  20 | -112.66675214
  21 | -112.66683834
  22 | -112.66688145
  23 | -112.66690301
  24 | -112.66691380
  25 | -112.66691919
-------------------------------------------------
Stopping because energy convergence was achieved.

-----+--------------
  #  | Energy
-----+--------------
   1 | -20.68609440
   2 | -11.39337409
   3 |  -1.56773075
   4 |  -0.79752832
   5 |  -0.64622306
   6 |  -0.64622306
   7 |  -0.54987945
   8 |   0.14682908
   9 |   0.14682908
  10 |   0.40542661

Outputting density files:
-------------------------
Pre-caching wavefunction values...[2.67012 seconds]
Writing 0001.dat
Writing 0002.dat
Writing 0003.dat
Writing 0004.dat
Writing 0005.dat
Writing 0006.dat
Writing 0007.dat
Writing 0008.dat
--------------------------------------------------------------
Done in 12.6541 seconds.
```
