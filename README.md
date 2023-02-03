# 3D genome reconstruction from partially phased Hi-C data
This repository contains files for the manuscript "3D genome reconstruction from partially phased Hi-C data".
<arxiv number>

## File descriptions

The repository contains three main parts:

* A Macualay2 script `unique_identifiability_conjecture.m2`, which contains the function SevenEq(nTrails) used test the unique identifiability conjecture (Conjecture 3.4 in the paper).

* A Maple script `ambig_identifiability_proof.mpl`, which does the rank computation referred to in the proof of identifiability in the fully ambiguous setting (Theorem 3.6 in the paper).

* A directory `reconstructions`, which contains files needed for the parts of the paper that concern concrete reconstructions (Sections 4 and 5). More specifically, it contains the following subdirectories:
  * `SNLC` containing all the MATLAB and Julia functions needed for the reconstruction method discussed in Section 4 of the paper. The most central functions are the Natlab function `estimate_disambiguated` (for estimating unambiguous loci), the Julia function `estimate_ambig_htpy` (for estimating ambiguous loci with homotopy continuation), the Matlab function `estimate_ambig`  (for refining estimations with local optimization), and unmix_chromosomes (for the clustering step). 
  * `synthetic_analysis`containing files needed for simulating and analyzing synthetic Hi-C data, including a Julia script `synthetic_example.jl` where these functions are used to produce Figures 3 (a)--(c) and S1 in the paper.
  * `patski_analysis` containing files needed for analyzing the patski dataset, including a Jupyter notebook `analysis_of_patski_data.ipynb` where the Figures 5, S2 and S3 are produced.
  
## Dependencies

The MATLAB function `estimate_disambiguated.m` relies on [ChromSDE version 2.2](http://glab.hzau.edu.cn/subpages/RESOURCES/softwares.php).
  
