% Feng et al "Polysaccharide utilization loci in Bacteroides determine population fitness and community-level interactions". 
% Submitted for publication in Cell Host & Microbe. Created by Yili Qian, Venturelli Lab, Nov 2020.
%%
clear all
close all
clc
warning off
%% Specify carbon source
carbon = 'pullulan';

%% MCMC setup
dyn_sys = @gLV_ParaVec;

logprior = @LogPrior;

num_chains = 1;

num_iti = 300000;   % number of samples generated

k_adapt = 0.2;  % adaptive stepping size. std of markov jump = k_adapt*std of previous jumps

burn_in = 10000;

%% set up experimental data
if strcmp(carbon,'pullulan') == 1
    exp_species = {1,2,3,4,5,[1 2],[1 3],[1 4],[1 5],6,[6 4],[6 5]};
    data_files = {'data_BU','data_AC','data_CC','data_ER','data_RI','data_BU_AC','data_BU_CC','data_BU_ER','data_BU_RI',...
        'data_BUD18','data_BUD18_ER','data_BUD18_RI'};
    cond_idx = {8,8,8,8,8,8,8,8,8,4,4,4};
    num_states = 6;
elseif strcmp(carbon,'glycogen') == 1
    exp_species = {1,2,3,4,5,[1 2],[1 3],[1 4],[1 5],6,[6 4]};
    data_files = {'data_BU','data_AC','data_CC','data_ER','data_RI','data_BU_AC','data_BU_CC','data_BU_ER','data_BU_RI',...
        'data_BUD18','data_BUD18_ER'};
    cond_idx = {2,2,2,2,2,2,2,2,2,1,1};
    num_states = 6;
elseif strcmp(carbon,'pectin') == 1
    exp_species = {1,2,3,4,5,[1 2],[1 3],[1 4],[1 5],6,[6 5]};
    data_files = {'data_BU','data_AC','data_CC','data_ER','data_RI','data_BU_AC','data_BU_CC','data_BU_ER','data_BU_RI',...
        'data_BUD18','data_BUD18_RI'};
    cond_idx = {7,7,7,7,7,7,7,7,7,3,3};
elseif strcmp(carbon,'laminarin') == 1
    exp_species = {1,2,3,4,5,[1 2],[1 3],[1 4],[1 5],6,[6 2],[6 4],[6 4],[6 5]};
    data_files = {'data_BU','data_AC','data_CC','data_ER','data_RI','data_BU_AC','data_BU_CC','data_BU_ER','data_BU_RI',...
        'data_BUD21','data_BUD21_AC','data_BUD21_CC','data_BUD21_ER','data_BUD21_RI'};
    cond_idx = {5,5,5,5,5,5,5,5,5,1,1,1,1,1};
    num_states = 6;
elseif strcmp(carbon,'glucomannan') == 1
    exp_species = {1,2,3,4,5,[1 2],[1 3],[1 4],[1 5]};
    data_files = {'data_BU','data_AC','data_CC','data_ER','data_RI','data_BU_AC','data_BU_CC','data_BU_ER','data_BU_RI'};
    cond_idx = {3,3,3,3,3,3,3,3,3};
    num_states = 5;
elseif strcmp(carbon,'pectic galactan') == 1
    exp_species = {1,2,3,4,5,[1 2],[1 3],[1 4],[1 5]};
    data_files = {'data_BU','data_AC','data_CC','data_ER','data_RI','data_BU_AC','data_BU_CC','data_BU_ER','data_BU_RI'};
    cond_idx = {6,6,6,6,6,6,6,6,6};
    num_states = 5;
elseif strcmp(carbon,'xyloglucan') == 1
    exp_species = {1,2,3,4,5,[1 2],[1 3],[1 4],[1 5]};
    data_files = {'data_BU','data_AC','data_CC','data_ER','data_RI','data_BU_AC','data_BU_CC','data_BU_ER','data_BU_RI'};
    cond_idx = {9,9,9,9,9,9,9,9,9};
    logprior = @LogPrior_xyloglucan;
    num_states = 6;
elseif strcmp(carbon,'inulin') == 1
    exp_species = {1,2,3,4,5,[1 2],[1 3],[1 4],[1 5]};
    data_files = {'data_BU','data_AC','data_CC','data_ER','data_RI','data_BU_AC','data_BU_CC','data_BU_ER','data_BU_RI'};
    cond_idx = {4,4,4,4,4,4,4,4,4};
    num_states = 5;
elseif strcmp(carbon,'glucose') == 1
    exp_species = {1,2,3,4,5,[1 2],[1 3],[1 4],[1 5]};
    data_files = {'data_BU','data_AC','data_CC','data_ER','data_RI','data_BU_AC','data_BU_CC','data_BU_ER','data_BU_RI'};
    cond_idx = {1,1,1,1,1,1,1,1,1};
    num_states = 5;
end

display(['carbon=' carbon ', num_iti=' num2str(num_iti) ', prior=' char(logprior)]);

%% initialization parameters
% initial step size for the burn-in period. Use adaptive stepping afterwards
delta_r = 0.02*ones(1,num_states);
delta_A = 0.02*ones(num_states,num_states);

step_min = 0.0005;

step_ini = [delta_r reshape(delta_A',1,num_states^2)];

% initial parameters
DetFitFile = ['para_seed_' carbon];

load(['MCMC_seed/' DetFitFile]);

%%
for q = 1:1:num_chains
    [para_vec,acceptance_vec,LogPosterior_vec] = AdaptiveMCMC_DynPara_v2(para0,logprior,dyn_sys,exp_species,...
        data_files,cond_idx,num_iti,burn_in,step_ini,k_adapt,step_min);

    acceptance_rate = movmean(acceptance_vec,50);     
    figure()
    subplot(1,2,1)
    plot(acceptance_rate,'k','linewidth',2)
    title('acceptance rate')
    
    subplot(1,2,2)
    plot(LogPosterior_vec,'k','linewidth',2)
    title('log likelihood')
     
    num_bars = 10;
    
    r_vec = para_vec(:,1:num_states);
    Aij_vec = para_vec(:,num_states+1:end);
    
    %% save
    FileName = ['MCMC_carbon_' carbon '_' datestr(now,'dd-mm-yyyy HH:MM:SS')];
    
    save(FileName,'para_vec','acceptance_vec','LogPosterior_vec','exp_species','num_iti','burn_in','step_ini',...
        'step_ini','acceptance_rate','num_states','r_vec','Aij_vec','acceptance_rate','k_adapt','cond_idx','data_files','exp_species');
end