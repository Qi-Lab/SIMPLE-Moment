# SIMPLE-Monent

This is an implementation of the SIMPLE-Moment method for ortho-positronium lifetime image reconstruction (paper title "High-resolution positronium lifetime tomography by the method of moments")

## Overview

- ./code: code for implementing the method of moments and results evaluation
- ./reconimg: reconstructed activity and lifetime image
- ./lm: list-mode data and time factors

## Usage

### Lifetime calculation

Run cal_lifetime_simple_moment.m for the SIMPLE-Moment method 

### Quantification and visualization of reconstructed lifetime images

Run eval_mean_and_var_mask.m and eval_mean_and_var_mask_moment.m to check the mean and s.d. of lifetime values and moment values respectively

### List-mode data

List-mode data with each event in int16 array [transaixal ID1, axial ID1, transaxial ID2, axial ID2, TOF] in folder ./lm/true. The correspoinding 1st to 3rd time weighting factors are in float32 arrays in folders ./lm/m1, ./lm/m2, and ./lm/m3

## Contact

For any questions or suggestions, feel free to contact Bangyan Huang at bybhuang@ucdavis.edu or Jinyi Qi at qi@ucdavis.edu
