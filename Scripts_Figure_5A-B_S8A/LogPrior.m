% Feng et al "Polysaccharide utilization loci in Bacteroides determine population fitness and community-level interactions". 
% Submitted for publication in Cell Host & Microbe. Created by Yili Qian, Venturelli Lab, Nov 2020.

function y = LogPrior(para)
% this function specifies the parameter priors

% find the number of species
N = (-1+sqrt(1+4*length(para)))*0.5;

% first N values are growth rates
R = para(1:N);

% rest are interaction coefficients -> matrix form
A_vec = para(N+1:end);
A = reshape(A_vec',N,N)';

% uniform prior
r_min = zeros(1,N);
r_max = 2*ones(1,N);

% prior of aii = unifpdf(a_min,0), prior of aij = unifpdf(a_min,a_max);
A_min = -2.5*ones(N,N);
A_max = 2.5*ones(N,N);

for i = 1:1:N
    A_max(i,i) = 0;
end

y = sum([log(unifpdf(R,r_min,r_max))';log(unifpdf(A(:),A_min(:),A_max(:)))]);
