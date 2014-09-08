3dmorse
==============

Oleg Devinyak (a,*), Dmytro Havrylyuk (b) and Roman Lesyk (b)

a)	Department of Pharmaceutical Disciplines, Uzhgorod National University, Uzhgorod 88000, Ukraine
	e-mail: o.devinyak@gmail.com
b)	Department of Pharmaceutical, Organic and Bioorganic Chemistry, Danylo Halytsky Lviv National Medical
	University, Lviv 79010, Ukraine

http://github.com/devinyak/3dmorse
--------------

3dmorse is an open-source small program to calculate 3D-MoRSE molecular descriptors. Currently it supports only MOPAC2012 output files (*.out) as input. The descriptors produced are 3D-MoRSE weighted with atomic mass, van der Waals volume, electronegativity, polarizability, atomic partial charge and unweighted descriptors. The naming convention is consistent with DRAGON 6 (powerful but only commercially available program for molecular descriptors calculation). That is, the numeration of 3D-MoRSE descriptors starts from 1, so, for example, Mor01u denotes unweighted descriptor with scattering parameter s=0 (since scattering parameter starts from zero), Mor02u denotes descriptor with  s=1 and so on. There is a possibility to obtain a table of 3D-MoRSE terms that correspond to each atomic pair in the molecular structure (for all descriptors at once). This table makes interpretation of 3D-MoRSE descriptors in a QSAR model much easier.

The typical usage of program is:
3dmorse path_to_input_file path_to_output_file <terms flag>
Terms flag is an optional argument, valid values are 0 (do not return 3D-MoRSE terms) or 1 (return 3D-MoRSE terms). The terms are not returned by default.

The output is comma separated values file with descriptors in columns. The output of terms has additional fragment "terms" in output file name and its columns are: N - serial number, firstAtom and secondAtom - correspond to atomic pair, s - scattering parameter, Distance - interatomic distance, term - corresponding summand value, weight - weighting scheme.

The program and its source are distributed under the GNU GPLv3 license.

For reference or citation use
Devinyak, O.; Havrylyuk, D.; Lesyk, R. 3D-MoRSE descriptors explained. Submitted to J. Mol. Graph. Model., 2014.
