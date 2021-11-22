% Feng et al "Polysaccharide utilization loci in Bacteroides determine population fitness and community-level interactions". 
% Submitted for publication in Cell Host & Microbe. Created by Yili Qian, Venturelli Lab, Nov 2020.

function [T,Y_mean,Y_std] = LoadFileFcn(dataFileName)
load(['data/' dataFileName])
T = time;
Y_mean = data.abs_mean;
Y_std = data.abs_std;