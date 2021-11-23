% Feng et al "Polysaccharide utilization loci in Bacteroides determine population fitness and community-level interactions". 
% Submitted for publication in Cell Host & Microbe. Created by Yili Qian, Venturelli Lab, Nov 2020.
clear all
close all
clc
%%
data_files = {'data_BU','data_AC','data_CC','data_ER','data_RI',...
    'data_BU_AC','data_BU_CC','data_BU_ER','data_BU_RI',...
    'data_BUD18_ER','data_BUD18_RI','data_BUD21_AC','data_BUD21_CC',...
    'data_BUD21_ER','data_BUD21_RI'};

species = {'BU','AC','CC','ER','RI','BUD18','BUD21'};

c_labels = [219 118 39  % AC
        42 157 72 % CC
        204 204 0   % ER
        255 10 10]; % RI
c_labels = c_labels/255;

BPB_mean_matrix = [];
BPB_std_matrix = [];
BPB_mean_sim_matrix = [];
butyrate_mean_vec = [];
butyrate_std_vec = [];
cond_vec = cell(1,2);

counter = 0;
for i = 1:length(data_files)
    FileName = data_files{i};
    load(['data/' FileName]);
    
    SpeciesName = FileName(6:end);
    
    Species_Idx = find(ismember(species,pairs));
    
    if max(Species_Idx)>5
        s = 'p';    % mutant
    else
        s = 's';    % WT
    end
    
    % find index of butyrate producers
    BPB_Idx = Species_Idx((Species_Idx>1)&(Species_Idx<6));
    
    num_fibers = length(data.carbon);
    
    for j = 1:num_fibers
        
        carbon = data.carbon{j};
        
        % load simulated abundance data. the simulated final abundance data 
        % "X_final" is generated using medium of the mcmc parameter outputs
        load(['MCMC_fit/para_' carbon])
        for q = 1:14
            II = find(X_final(q,:) ~= 0 & X_final(q,:) ~= 0);
            if II == Species_Idx
                row_idx = q;
            end
        end
        
        butyrate_mean = data.met_mean{j};
        butyrate_std = data.met_std{j};
        
        BPB_mean = data.abs_mean{j}(end,[2:5]);
        BPB_std = data.abs_std{j}(end,[2:5]);

        BPB_mean_sim = X_final(row_idx,[2:5]);
        
        if isnan(butyrate_mean)~=1
            counter = counter+1;
            
            cond_name{counter,1} = [FileName carbon];
            
            butyrate_mean_vec = [butyrate_mean_vec;butyrate_mean];
            butyrate_std_vec = [butyrate_std_vec;butyrate_std];
            
            BPB_mean_matrix = [BPB_mean_matrix;BPB_mean];
            BPB_std_matrix = [BPB_std_matrix;BPB_std];
            
            BPB_mean_sim_matrix = [BPB_mean_sim_matrix;BPB_mean_sim];
            
            cond_vec{counter,1} = pairs;
            cond_vec{counter,2} = data.carbon{j};
            
            c_vec(counter,:) = c_labels(BPB_Idx-1,:);
            
            s_vec{counter} = s;
            
        end
    end
end

% fit a linear function butyrate = \sum_i k_i * butyrate producer i abundance
LM = fitlm(BPB_mean_matrix,butyrate_mean_vec,'Intercept',false);
k_model = LM.Coefficients.Estimate;
butyrate_predicted_vec = BPB_mean_matrix*k_model;
butyrate_predicted_std_vec = BPB_std_matrix*k_model;

figure(2)
for i = 1:length(butyrate_predicted_vec)
    if length(cond_vec{i,1}) == 2   % coculture
        scatter(butyrate_predicted_vec(i),butyrate_mean_vec(i),50,'k',s_vec{i},'linewidth',1,'markerfacecolor',c_vec(i,:));
    else  % monoculture
        scatter(butyrate_predicted_vec(i),butyrate_mean_vec(i),50,c_vec(i,:),s_vec{i},'linewidth',1);
    end
    hold on
    errorbar(butyrate_predicted_vec(i),butyrate_mean_vec(i),-butyrate_std_vec(i),butyrate_std_vec(i),...
        -butyrate_predicted_std_vec(i),butyrate_predicted_std_vec(i),...
        'color',c_vec(i,:),'linewidth',0.2)
    hold on
end

plot([0 35],[0 35],'--k','linewidth',1)

xlabel('predicted butyrate (mM)')
ylabel('actual butyrate (mM)')
print(gcf,'Butyrate.pdf','-dpdf','-bestfit','-painters')
set(gca,'fontsize',10)
xlim([0 35]);ylim([0 35])

set(gcf,'position',[0 0 300 300])