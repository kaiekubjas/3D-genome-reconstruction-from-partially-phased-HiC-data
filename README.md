# Repository for the manuscript "3D genome reconstruction from partially phased Hi-C data"
This repository contains files for the manuscript "3D genome reconstruction from partially phased Hi-C data".
<arxiv number>

## File descriptions

The repository contains three main parts:

* A Macualay2 script `unique_identifiability_conjecture.m2`, which contains the function SevenEq(nTrails) used test the unique identifiability conjecture (Conjecture 3.4 in the paper).

* A Maple script `ambig_identifiability_proof.mpl`, which does the rank computation referred to in the proof of identifiability in the fully ambiguous setting (Theorem 3.6 in the paper).

* A directory `reconstructions`, which contains files needed for the parts of the paper that concern concrete reconstructions (Sections 4 and 5). More specifically, it contains the following subdirectories:
  * `ChromSDE_program2.2` containing files from the software package ChromSDE.
  * `SNLC` containing all the MATLAB and Julia functions needed for the reconstruction method discussed in Section 4 of the paper.
  * `synthetic_analysis`containing files needed for simulating and analyzing synthetic Hi-C data, including a Julia script `synthetic_example.jl` where these functions are used to produce Figures 3 (a)--(c) and S1 in the paper.
  * `patski_analysis` containing files needed for analyzing the patski dataset, including a Jupyter notebook `analysis_of_patski_data.ipynb` where the Figures 5, S2 and S3 are produced.
  
