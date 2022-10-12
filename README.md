# Please see LICENSE.md as well as the journals of the publications below for copyright information.

Description of codes:

1. SYK_corr_func.nb

This Mathematica notebook studies a model of condensed-matter called the SYK model. This model describes the quantum mechanics of a Fermi particle with a 4-point interaction. See https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.106002 for a review. We study the Hamiltonian of this model. It is a product of 4 matrices obeying certain constraints. They are constructed using Array functions in Mathematica, the Hamiltonian is then calculated and diagonalised. The resulting eigenvalues and eigenvectors are used to calculate two-point function and four-point function of the Fermions. These quantities are important for the physics of SYK model. This code was used in the publication: Jan de Boer, Rik van Breukelen, Sagar F. Lokhande, Kyriakos Papadodimas, Erik Verlinde; Probing Typical Black Hole Microstates; https://doi.org/10.1007/JHEP01%282020%29062 JHEP 01 (2020) 062.

2. SYK_matlab.m

MatLab notebook that once again studies SYK model. It basically follows similar logic as in SYK_corr_func.nb, but with matrices of higher dimensions. The Fermion matrices are defined, Hamiltonian calculated and diagonalised, then the eigenvalues and eigenvectors used to calculate correlation functions. 

3. random_matrices_tfd.nb

This Mathematica code studies random matrices, in particular the distribution of their eigenvalues as the dimension of the matrix grow very large. The model contains two Hamiltonians (large-dimensional symmetric real matrices) along with an interaction matrix which is much more general. We study specific cases for the interaction matrix. Of particular importance in Physics is the ground state of this model, called the Thermofield Double State. The code first constructs the component matrices and then the model, bootstraps it (diagonalizing the finding eigenvalues and eigenvectors) and finally analyses the resulting data. There are many scatter plots, histograms, density graphs and other fun stuff relevant in the study of random matrices. This notebook was used in the following publication:
William Cottrell, Ben Freivogel, Diego M. Hofman, Sagar F. Lokhande; How to Build the Thermofield Double State; https://doi.org/10.1007/JHEP02%282019%29058 JHEP 02 (2019) 058.

4. my_ls.py

This is a command line Python script which defines a Python alternative to the Unix Shell "ls" command. The "ls" command lists the contents of a directory, and comes with many useful options. This Python script uses the glob package to mimic the same functionality.
