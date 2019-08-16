# Representational Similarity Analysis tutorial

Author: Choong-Wan Woo (Sungkyunkwan University) https://cocoanlab.github.io/
Date: 2019/8/17 @ KHBM 2019 Summer school

## Slides

Download: You can download the slide PDF [here](/slide/Woo_KHBM2019_summerschool_RSA.pdf)

## Dependencies

To run the **Matlab scripts** `tutorial_main.mlx`, or `tutorial_main.m`, you will need the following tools installed in your computer. The code and results can be viewed in `tutorial_main.html` or `tutorial_main.pdf`. 

+ Matlab (> 2016 version)

+ [CanlabCore tools](https://github.com/canlab/CanlabCore)

+ For full functionality, make sure to install:
	- Matlab Statistics and Machine Learning toolbox
	- [Statistical Parametric Mapping (SPM) toolbox](https://www.fil.ion.ucl.ac.uk/spm/) 

The Matlab script `tutorial_main.mlx` was tested on macOS High Sierra with Matlab 9.5 R2018b.


See also RSA toolbox written by Kriegeskorte group: [https://github.com/rsagroup/rsatoolbox](https://github.com/rsagroup/rsatoolbox)

## Dataset
from Woo et al., 2014, Nat Comms; Paper PDF can be downloaded [here](https://cocoanlab.github.io/pdfs/Woo_2014_NatComms.pdf)

Participants: _N_ = 59

#### Tasks
There were two types of tasks, and in each task, there were two conditions (2 x 2 design)
- Physical pain task (Heat, Warmth conditions)
- Social pain task (Rejection, Friends conditions)

## Analysis plan
- Step 1: Computing RDMs for each participant, for each region (4 ROIs: aINS, dACC, S2/dpINS, TPJ),  and visualize the RDMs
- Step 2: Comparing brain and model RDMs: We have two model RDMs: physical vs. social, aversive vs. non-aversive, and 4 ROIs
- Step 3: Statistical inference: Using the neurosync mask that Woo 2014 used, we will do statistical inference
