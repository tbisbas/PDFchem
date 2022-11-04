The PDFchem algorithm is a continuation of the work by Bisbas et al. (2019, MNRAS, 485, 3097). PDFchem is a combined Fortran 90/95 and Python algorithm. In principle, the Python script (written for Jupyter notebook) is a plotting tool in which PDFchem can also be executed.

To run PDFchem you need first to download the PDR simulations [LINK PENDING]

Once downloaded, you will need to untar it using the command
tar xvzf simulations.tgz
This will extract four directories in which all PDR files will be placed. The PDR files are in ASCII format for an easy use in case of further/other studies the user wishes to do. For more information about the files and how they are strutured, please refer to [THIS ARTICLE].

Next, you will need to compile the pdfchem_algorithm.f90 program. To do this, you need to have fortran compilers and type:
gfortran -o pdfchem_algorithm pdfchem_algorithm.f90

Once the algorithm is compiled, all you need to do is to load the PDFchem.ipynb file using:
jupyter notebook PDFchem.ipynb
The notebook is self-explanatory. All calculations are performed for a log-normal distribution only. For any other distributions you will need to do the relevant edits and modifications in the pdfchem_algorithm.f90 file. In the notebook you can insert your own mean (Av) and the widht (sigma) of the Av-PDF you wish to consider. 

PDFchem will then perform calculations and it will plot some default figures. If you wish to plot something more specific, you will need to edit the relevant parts using the variable of the species/emittor available. Currently, the available abundances are (copied from the PDFchem.ipynb file):

#Sequence of species
# 1,H3+;  2,He+;   3,Mg;     4,H2+;   5,CH5+
# 6,O2;   7,CH4+;  8,O+;     9,OH+;  10,Mg+
#11,C+;  12,CH4;  13,H2O+;  14,H3O+; 15,CO+
#21,CH;  22,CH3;  23,HCO+;  24,CH2+; 25,C
#26,He;  27,CH+;  28,CO;    29,OH;   30,O
#31,H2;  32,H;    33,e-; 

The available brightness temperatures consider the carbon cycle only:
[CII] 158um
[CI] (1-0), (2-1)
CO (1-0) ... (10-9)

For any question, please do not hesitate to contact me at tbisbas@gmail.com

Thomas Bisbas
thomasbisbas.com
