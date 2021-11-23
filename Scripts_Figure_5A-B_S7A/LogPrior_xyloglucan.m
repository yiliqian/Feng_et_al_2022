function y = LogPrior_xyloglucan(para)

% find the number of species
N = (-1+sqrt(1+4*length(para)))*0.5;

% first N values are growth rates
R = para(1:N);

% rest are interaction coefficients -> matrix form
A_vec = para(N+1:end);
A = reshape(A_vec',N,N)';

% normal prior
R0 = [0.14	0.07	0.003	0.017	0.12	0];
A0 = [-0.05	0	0	0	-0.3	0
    0.06	-0.22	0	0	0	0
    0.18	0	-0.9	0	0	0
    0.012	0	0	-0.24	0	0
    0.11	0	0	0	-0.77	0
    0	0	0	0	0	0];

std_R = max(R0*0.1,0.05);
std_A = max(abs(A0)*0.1,0.05);

r_min = zeros(1,N);
r_max = 2*ones(1,N);

y = sum([log(normpdf(R,R0,std_R))';log(normpdf(A(:),A0(:),std_A(:)))])+sum(log(unifpdf(R,r_min,r_max)));
