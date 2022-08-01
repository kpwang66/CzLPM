# init

## Description

An algebraic model is devised to predict oxygen concentration in molten silicon (produced using a method known as Czochralski method) given a number of physical parameters about the system.  The model is built on physical arguments, though three parameters are difficult to assess experimentally.  Thus, I used the optimization toolbox to fit the model to a few limiting cases in the experimental data to deduce the experimental values of these three parameters.

For an in-depth understanding of model derivation, please see the included PDF or visit here: [https://doi.org/10.1016/j.jcrysgro.2021.126384](https://doi.org/10.1016/j.jcrysgro.2021.126384)

### Abstract

Lumped-parameter models are derived from boundary layer and other physical arguments to describe oxygen concentration levels during the Czochralski (CZ) growth of silicon. These models are assessed against predictions from a detailed, high-fidelity 2D-3D numerical simulation of the entire CZ puller, whose solutions are realistic but require intense computational effort. Comparisons of predictions show that the lumped-parameter model captures the correct trends of melt oxygen levels influenced by melt height, crucible rotation, and crystal rotation. A simple fitting of coefficients provides reasonably good quantitative predictions by the lumped-parameter model, and its near-instantaneous computations make it an interesting candidate for real-time growth optimization and control. Possible model improvements and extensions are discussed.

### Relevant files:
- `CzO_Data.mat` - MATLAB data file containing data from numerous Czochralski direct numerical simnulations
- `CzROM_newmono_v5.m` - Main MATLAB file to run
- `Physically-based, lumped-parameter models for the prediction of oxygen concentration during Czochralski growth of silicon crystals.pdf`

## Setup and Use

Make sure CzO_Data.mat is present and in the MATLAB path.  The entire model is executed from the `CzROM_newmono_v5.m` MATLAB script.  The script reproduces the calculations and plots found in Figures 2 and 3 of the above publication.