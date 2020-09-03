# patlm

Peptide adsorption to lipid membranes: computational model

## Description
This is the source code of the computational model used in:

Pedro G. Ramírez, Mario G. Del Popolo, Jorge A. Vila and Gabriel S. Longo "Thermodynamics of Cell Penetrating Peptides on Lipid Membranes: Sequence and Membrane Acidity Regulate Surface Binding" (under review)

## Instructions

You will need a Fortran 95 compiler and a functioning installation of [KINSOL](https://computing.llnl.gov/projects/sundials/kinsol). To compile use ```compi.sh```. Input the path to KINSOL in that file.

You can try the example by running the compiled executable from the example folder. All required input files are contained inside the example folder.
