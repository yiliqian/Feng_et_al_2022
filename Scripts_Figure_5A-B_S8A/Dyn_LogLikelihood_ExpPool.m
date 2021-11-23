% Feng et al "Polysaccharide utilization loci in Bacteroides determine population fitness and community-level interactions". 
% Submitted for publication in Cell Host & Microbe. Created by Yili Qian, Venturelli Lab, Nov 2020.

function LogLikelihood = Dyn_LogLikelihood_ExpPool(exp_species,data_files,cond_idx,para_current,dyn_sys)
% Dyn_LogLikelihood_ExpPool(exp_species,data_files,cond_idx,para_current,dyn_sys)
% This script pools together the likelihood of a set of parameter over all experimental
% conditions. Let L be the number of experimental data sets.
% Input:
% - exp_species: 1 x L cell. Each cell contains a vector 1 x num_species in
% size that specifies the columns in each data file (which corresponds to
% species) that should be read.
% - data files: 1 x L cell. Each cell contains a string specifies the
% experimental data files to read.
% - cond_idx: 1 x L cell. Which conditions needs to be read in each data
% file. Each data file should contain a cell (1xnum_conditions size). Each
% cell contains the data for that specific doncition.
% - para_current: 1xp row vector, sampled parameter
% - dyn_sys: dynamical system to simulation, a function handle

% figure out the number of states
N = (-1+sqrt(1+4*length(para_current)))*0.5;

% number of experimental data sets to pool together
L = length(exp_species);

% put parameters in matrix format
r = para_current(1:N);
A = reshape(para_current(N+1:end),N,N)';

% cumulative likelihood
LogLikelihood_cum = 0;

for i = 1:1:L
    species_idx = exp_species{i};
    dataFileName = data_files{i};
    
    % parameters needed for this experiment
    r_exp = r(species_idx);
    A_exp = A(species_idx,species_idx);
    
    [nn,~] = size(A_exp);
    
    A_exp_vec = reshape(A_exp',1,nn^2);
    
    para_exp_vec = [r_exp,A_exp_vec];
    
    % how many experimental conditions needs to be used
    num_conditions = length(cond_idx{i});
    
    % this loads the data from all experimental conditions
    [T,Y_mean,Y_std] = LoadFileFcn(dataFileName);
    
    for j = 1:1:num_conditions
        
        y_mean = Y_mean{cond_idx{i}(j)}(:,species_idx);
        y_mean = max(y_mean,0,'includenan');
        y_std = Y_std{cond_idx{i}(j)}(:,species_idx);
        y_std = max(y_std,0.02);

        % use the first measurement as initial condition
        X0 = y_mean(1,:);
        
        LogLikelihood_cum = LogLikelihood_cum + Dyn_LogLikelihood_SingleExp(para_exp_vec,dyn_sys,T,species_idx,X0,y_mean,y_std);
    end
end
LogLikelihood = LogLikelihood_cum;