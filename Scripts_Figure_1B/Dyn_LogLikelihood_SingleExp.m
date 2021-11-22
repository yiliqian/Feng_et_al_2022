% Feng et al "Polysaccharide utilization loci in Bacteroides determine population fitness and community-level interactions". 
% Submitted for publication in Cell Host & Microbe. Created by Yili Qian, Venturelli Lab, Nov 2020.

function LogLikelihood = Dyn_LogLikelihood_SingleExp(para_current,dyn_sys,time,X0,Y_mean,Y_std)
% Dyn_LogLikelihood_SingleExp(para_current,dyn_sys,time,X0,Y_mean,Y_std)
% compute likelihood of a sample parameter para_current given measurement Y,
% measurement noise Y_std, sampling time
% Input:
% - para_current: parameter sampled [1 x p vector]
% - dyn_sys: simulated function
% - time: sampling time [N x 1 vector]
% - species_idx: take into account the dynamics of which species specified
% by vector
% - X0: initial condition [n x 1 vector]
% - Y_mean: experimentally meausred mean [N x num_outputs]
% - Y_std: experimentally measured std [N x num_outputs]

L = length(time);

[~,X] = ode23s(@(t,x) dyn_sys(t,x,para_current),time,X0);

% % discard states not included in species_idx (e.g., when these species are not included in experiments,
% % so mean and std are both 0).
% X = X(:,species_idx);

% if system unstable set trajectory to infity
[L_sim,W_sim] = size(X);

if L_sim < L
    X = inf*ones(L,W_sim);
end

p_y_theta = zeros(L-1,1);   % compute likelihood of each temporal data point, remove initial point since it is deterministic

for i = 1:1:L-1
    for j = 1:1:W_sim
        % evaluate likelihood: assume channels are independent
        p_y_theta(i) = p_y_theta(i)+lognormalpdf(X(i+1,j),Y_mean(i+1,j),Y_std(i+1,j));
    end
end

LogLikelihood = sum(p_y_theta);

function y = lognormalpdf(x,y_mean,y_std)

if sum(isnan(y_mean)) == 0
    y = log(normpdf(x,y_mean,y_std));
else
    y = 0;
end