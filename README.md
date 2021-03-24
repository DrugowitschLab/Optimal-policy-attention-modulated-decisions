# Optimal-policy-attention-modulated-decisions


This repository contains the raw data and MATLAB scripts to reproduce the main and supplemental figures of the following article:

Jang, A. I., Sharma, R., & Drugowitsch, J. (2021). Optimal policy for attention-modulated decisions explains human fixation behavior. eLife (In press).

## Installation
To run the scripts, download & extract them into a folder of your choice, and navigate to this folder within MATLAB.

The 'data' folder contains the human behavioral data from the following two papers (shared here with authors' permission):
1. datastruct_aVBDM.mat: Krajbich, I., Armel, C., & Rangel, A. (2010). Visual fixations and the computation and comparison of value in simple choice. Nature neuroscience, 13(10), 1292-1298.
2. datastruct_aPDM.mat: Tavares, G., Perona, P., & Rangel, A. (2017). The attentional drift diffusion model of simple perceptual decision-making. Frontiers in neuroscience, 11, 468.

The 'data/optimal_policy_save' folder contains pre-computed optimal policy spaces based on different parameter values. If these policy spaces are deleted, they will be re-computed and saved when you run the script. 

The 'functions' folder contains all the custom MATLAB functions that are used in the main scripts. This directory will be added in the beginning of each script. 

All scripts have been tested and run under MATLAB R2020a on MacOS Catalina. They should also work on other operating systems or later MATLAB versions. There are two scripts, which should utilize the functions included in this package to plot all the figures. Simply set the 'rootdir' variable to the directory that includes the scripts, and they should run correctly.
1. Jang2021_Figures.m
2. Jang2021_Figures_supplement.m

## Abstract
Traditional accumulation-to-bound decision-making models assume that all choice options are processed with equal attention. In real life decisions, however, humans alternate their visual fixation between individual items to efficiently gather relevant information (Yang et al., 13 2016). These fixations also causally affect one’s choices, biasing them toward the longer-fixated item (Krajbich et al., 2010). We derive a normative decision-making model in which attention enhances the reliability of information, consistent with neurophysiological findings (Cohen and Maunsell, 2009). Furthermore, our model actively controls fixation changes to optimize information gathering. We show that the optimal model reproduces fixation-related choice biases seen in humans and provides a Bayesian computational rationale for this phenomenon. This insight led to additional predictions that we could confirm in human data. Finally, by varying the relative cognitive advantage conferred by attention, we show that decision performance is benefited by a balanced spread of resources between the attended and unattended items.
