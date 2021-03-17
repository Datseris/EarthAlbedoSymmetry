# Earth's albedo and its symmetry: code base.

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project accompanying our paper *Earth's albedo and its symmetry*, authored by George Datseris and Bjorn Stevens.

To locally reproduce this project, do the following:

0. Download this code base as is.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```
   This will install all necessary packages for you to be able to run the scripts and everything should work out of the box.
2. This project uses CERES data, which are cited properly within our paper. Specifically, we use monthly averages of the EBAF data product, both TOA and SFC versions, as well as the ice + ocean area fraction from the SYN1deg datasets. These datasets are present in the `data` folder.
3. All code that makes all figures of our paper is in the Jupyter notebook file `figures.ipynb`. This file is guaranteed`*` to run out of the box. You can also see the file using nbviewer here: https://nbviewer.jupyter.org/github/Datseris/EarthAlbedoSymmetry/blob/main/figures.ipynb

`*` If the Jupyter notebook or the Julia script `papers/figures.jl` do not run out of the box, please open an Issue. Should be an easy fix.
