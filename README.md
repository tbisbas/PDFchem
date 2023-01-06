# PDFchem v2.0 (2022)

This version is based on the paper by [Bisbas et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.519..729B/abstract). It is a continuation of the work by [Bisbas et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.3097B/abstract). PDFchem is a combined Fortran 90/95 and Python algorithm. In principle, the Python script (written for Jupyter notebook) is a plotting tool in which PDFchem can also be executed and it serves as a wrapper.

If you use PDFchem for your research, please cite the above papers and [the code](https://ui.adsabs.harvard.edu/abs/2022ascl.soft11014B/abstract).

## Installing PDFchem

To run PDFchem you need first to download the PDR simulations from [this Zenodo link](https://zenodo.org/record/7310833).

Once downloaded, you will need to untar it using the command
```
tar xvzf simulations.tgz
```
This will extract four directories in which all PDR files will be placed. The PDR files are in ASCII format for an easy use in case of further/other studies the user wishes to do. For more information about the files and how they are strutured, please refer to the [3D-PDR manual](https://uclchem.github.io/3DPDR_manual.pdf).

Next, you will need to compile the pdfchem_algorithm.f90 program. This will be done only once. To do this, you need to have fortran compilers and type:
```
gfortran -o pdfchem_algorithm pdfchem_algorithm.f90
```

## Running PDFchem

Once the algorithm is compiled, all you need to do is to load the PDFchem.ipynb file using:
```
jupyter notebook PDFchem.ipynb
```
The notebook is self-explanatory. All calculations are performed for a log-normal distribution only. For any other distributions you will need to do the relevant edits and modifications in the pdfchem_algorithm.f90 file. In the notebook you can insert your own mean (Av) and the width (sigma) of the Av-PDF you wish to consider. 

PDFchem will then perform calculations and it will plot some default figures. If you wish to plot something more specific, you will need to edit the relevant parts using the variable of the species/emittor available. The list of available species can be found in the PDFchem.ipynb file.

The available brightness temperatures consider the carbon cycle only:
- [CII] 158um
- [CI] (1-0), (2-1)
- CO (1-0) ... (10-9)

## Interactive plots

If you are interested to explore the results with interactive plots, you can use the following notebook (with thanks to Theodoros Topkaras for the contribution)
```
jupyter notebook PDFchem_plotly.ipynb
```

### Contact

For any question, please do not hesitate to contact me at tbisbas@gmail.com or tbisbas@zhejianglab.com

Thomas Bisbas ([thomasbisbas.com](http://thomasbisbas.com))
