% Feng et al "Polysaccharide utilization loci in Bacteroides determine population fitness and community-level interactions". 
% Submitted for publication in Cell Host & Microbe. Created by Yili Qian, Venturelli Lab, Nov 2020.

clear all
close all
clc
warning off
%%
carbon_vec = {'Arabinogalactan','Galactanmannan','Glucomannan',...
    'Glycogen','Inulin','Laminarin','PecticGalactan','Pectin',...
    'Pullulan','TypeIIMucin','Xyloglucan','Glucose'};

species_vec = {'BU','delPUL6','delPUL7','delPUL11','delPUL1143',...
    'delPUL12','delPUL13','delPUL16','delPUL17','delPUL18',...
    'delPUL21','delPUL2223','delPUL24','delPUL28','delPUL32',...
    'delPUL34','delPUL35','delPUL37','delPUL43','delPUL44',...
    'delPUL47','delPUL48','delPUL49','delPUL54'};

%% load experimental data
for QQ = 1:length(species_vec)
    species = species_vec{QQ};
    for q = 1:length(carbon_vec)
        carbon = carbon_vec{q};
        
        DataFileName = ['data/' species '_' carbon];
        
        load(DataFileName);
        
        [MAX,I_max] = max(MEAN);
        
        I_small = find(MEAN<= 0.8*MAX);
        I_rmv = I_small(find(I_small > I_max));
        if min(I_rmv)<=3
            time = time(1:3);   % ensure at least 3 data points are left (only 2 for fitting as first is IC)
            MEAN = MEAN(1:3);
            STD = STD(1:3);
        else
            time(I_rmv) = [];
            MEAN(I_rmv) = [];
            STD(I_rmv) = [];
        end
        
        STD = max(STD,0.01);
        %% load deterministic fitting data
        load('SeedPara')
        
        strain_idx = find(strcmp(species,strain_vec) == 1);
        carbon_idx = find(strcmp(carbon,carbon_vec) == 1);
        
        R_det = R_opt(strain_idx,carbon_idx);
        A_det = A_opt(strain_idx,carbon_idx);
        
        %% MCMC setup
        dyn_sys = @gLV_ParaVec;
        
        num_states = 1; % total number of species considered (i.e., set from which subcommunities are picked)
        
        num_iti = 100;   % number of samples generated
        
        k_adapt = 0.05;  % adaptive stepping size. std of markov jump = k_adapt*std of previous jumps
        
        burn_in = 0.1*num_iti;
        
        display(['species=' species ', carbon=' carbon ', num_iti=' num2str(num_iti)]);
        %% initialization parameters
        % initial step size for the burn-in period. Use adaptive stepping afterwards
        delta_r = 0.02;
        delta_A = 0.02;
        
        step_min = 0.0005;
        
        step_ini = [delta_r delta_A];
        
        para0 = [R_det A_det];
        
        [para_vec,acceptance_vec,LogPosterior_vec] = AdaptiveMCMC_DynPara_v2(para0,dyn_sys,...
            time,MEAN,STD,num_iti,burn_in,step_ini,k_adapt,step_min);
        
        acceptance_rate = movmean(acceptance_vec,50);
        
        r_vec = para_vec(:,1:num_states);
        Aij_vec = para_vec(:,num_states+1:end);
        
        %% save
        FileName = ['MCMC_run/MCMC_' species '_' carbon];
        
        save(FileName,'para_vec','acceptance_vec','LogPosterior_vec','num_iti','burn_in','step_ini',...
            'step_ini','num_states','r_vec','Aij_vec','acceptance_rate','k_adapt','time',...
            'MEAN','STD');
    end
end