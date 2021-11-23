% Feng et al "Polysaccharide utilization loci in Bacteroides determine population fitness and community-level interactions". 
% Submitted for publication in Cell Host & Microbe. Created by Yili Qian, Venturelli Lab, Nov 2020.

function [para_vec,acceptance_vec,LogPosterior_vec] = AdaptiveMCMC_DynPara_v2(para0,prior,dyn_sys,exp_species,...
    data_files,cond_idx,num_iti,burn_in,step_ini,k_adapt,step_min)
% Dyn_LogLikelihood_ExpPool(exp_species,data_files,cond_idx,para_current,dyn_sys) 
% Adaptive MCMC to identify dynamical system parameters.
% Need to have the following files in the same directory:
% 1. dyn_sys.m [dynamical sys to simulate]
% 2. LogPrior.m [parameter prior distribution]
%%%%%%%%%%%%%
% Inputs:
% para0 - initial guess of the parameters [row vector 1 x p, p = # of parameters]
% dyn_sys - dynamical system to simulate [use function handle to specify]
% time - measurement time points [column vector N x 1, N = # of temporal points]
% X0 - initial condition of dynamical system [column vector 1xn, b = # of states]
% mean_meas - measured mean value [matrix N x n]
% std_meas - measured standard deviation [matrix N x n]
% num_iti - total number of iterations [positive integer]
% burn_in - # of iterations to discard as burn-in [poisitve integer << num_iti]
% step_ini - initial std of proposal (normal) distribution before
% adaptation begines [row vector 1xp]
%%%%%%%%%%%%
% Outputs:
% para_vec - (parameter) sample trajectory from MCMC [matrix num_iti x p]
% acceptance_vec - whether a (parameter) sample is accepted (1) or rejected (0) [vector num_itix1]
% LogPosterior_vec - loglikelihood of each sample
% X_vec - simulated trajectory for each sample

%% set up
% get the number of parameters
num_parameters = length(para0);

% minimum step size for adaptive MCMC
% step_min = 0.005;  

% records log likelihood (start with index 0 to always accept the first sample)
LogPosterior_vec = -inf*ones(num_iti+1,1);

% records if a sample is accepted or not
acceptance_vec = zeros(num_iti,1);

% acceptance probability
a_vec = zeros(num_iti+1,1);

% sampled r vector
para_vec = zeros(num_iti+1,num_parameters);

% generate uniform sample from [0,1] to compare with acceptance ratio take
% log because the likelihoods are also in log
seed = log((rand(num_iti,1)));

% generate normal sample from [0,1] to perform Markovian jump
para_seed_jump = randn(num_iti,num_parameters);

%% compute loglikelihood of para0
% loglikelihood associated with dynamic measurements
LogLikelihood_dyn = Dyn_LogLikelihood_ExpPool(exp_species,data_files,cond_idx,para0,dyn_sys);

% compute log posterior
LogPosterior = LogLikelihood_dyn + LogPrior(para0);

LogPosterior_vec(1) = LogPosterior;

%% run Markov chain
para_vec(1,:) = para0;

% current state of parameters
para_current = para0;

% set jump std during burn-in to be the initial step size
delta_para = step_ini;

for j = 2:1:num_iti+1
    
    % display progress every 100 jumps
    if mod(j,100) == 0
        display(['iteration ' num2str(j)])
    end
    
    % adaptive MCMC scheme from Robserts and Rosenthal 2008
    if j>burn_in
%         para_std = 2.38^2*std(para_vec(1:j-1,:))/num_parameters;
        para_std = k_adapt*std(para_vec(1:j-1,:))/num_parameters;
        
        % use standard deviation from existing samples to set new "jump size"
        delta_para = max(para_std,step_min);
    end
    
    % propose next jump state
    para_current = para_current + para_seed_jump(j-1,:).*delta_para;
    
    % compute log likelihood
    LogLikelihood_dyn = Dyn_LogLikelihood_ExpPool(exp_species,data_files,cond_idx,para_current,dyn_sys);
    
    % compute log posterior
    LogPosterior = LogLikelihood_dyn + prior(para_current);
    
    % log acceptance ratio
    a = min(0,LogPosterior-LogPosterior_vec(j-1));
    a_vec(j-1) = a;
    
    if a > seed(j-1)
        % if greater than seed, accept, record parameter sample and likelihood
        acceptance_vec(j-1) = 1;
        
        para_vec(j,:) = para_current;

        LogPosterior_vec(j) = LogPosterior;
    else
        % if less than seed, reject, do not record parameter, replace
        % likelihood with previous likelihood
        acceptance_vec(j-1) = 0;
        
        para_vec(j,:) = para_vec(j-1,:);
        
        para_current = para_vec(j-1,:);
        
        LogPosterior_vec(j) = LogPosterior_vec(j-1);
    end
end