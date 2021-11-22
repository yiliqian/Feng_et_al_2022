README:

This folder contains scripts and all raw data used to generate Figure 1B in the manuscript: 

Feng et al "Polysaccharide utilization loci in Bacteroides determine population fitness and community-level interactions", submitted for publication in Cell Host & Microbe.

A brief description of each items in this folder can be found below:

- "MCMC_ini.m": MATLAB script to set up MCMC parameters (number of iterations, prior functions, adaptive gain parameters etc.) and initiate an MCMC run. To run parameter inference, use this script. 

- "data" folder: A folder containing monoculture experimental data for parameter inference. 

- "SeedPara.mat": A MATLAB data file that stores the initial conditions for MCMC runs used for the paper. Different initial conditions could be used. 

****************
 
MATLAB functions below are all subfunctions for "MCMC_ini.m". 

- "gLV_ParaVec.m":  MATLAB function to simulated the generalized Lotka-Volterra equation.

- "AdaptiveMCMC_DynPara_v2.m": MATLAB function to perform adaptive MCMC.

- "Dyn_LogLikelihood_ExpPool.m": MATLAB function to compute the likelihood of a proposed parameter given a pool of experimental data (i.e., different experimental data sets, including co-culture data and monoculture).

- "Dyn_LogLieklihood_SingleExp.m": MATLAB function to compute the likelihood of a proposed parameter given a single experimental data trajectory.

- "LogPrior.m": MATLAB function to specify the parameter priors for MCMC. This uniform prior is used for all carbon sources except xyloglucan.

- "LogPrior_xyloglucan.m": MATLAB function to specify the parameter prior for parameter inference in the xyloglucan data set. (See details in STAR method and Supplemental Data 6.)

- "LoadFileFcn.m": MATLAB function to load experimental data for parameter inference. 

  
