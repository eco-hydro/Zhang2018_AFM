clc
%% the number of free parameters 

% In their paper, acording to Fig.1, also known that the entire time series 
%   (2003-2007) (Page404, line7)was filterd in a single run using whittaker smoother,
%   so the the length of input data Y is (46 composites per year) * 5 year = 230; 
m = 46*5; 
d = 3;           % in their paper, third order difference was used (Eq. 13)
lambda = 2;      % lambda = 2 or 15

%% the first solution to calculate the number of free parameter
%   According to: Eilers, P.H.C., & Marx, B.D.(1996). Flexible smoothing 
%   with B-splines and penalties. Statistical Science, 11(2), 89–121.
E = speye(m); %sp
D = diff(E, d);
C = chol(E + lambda * D' * D);

Qb = C' * C;
Qlambda = lambda * (D' * D);
freedom = trace((Qb + Qlambda)\Qb)
% freedom = trace((E + Qlambda)\E);

%% the second
%   According to: Eilers, P. H. C. (2003). A perfect smoother. Analytical Chemistry, 75(14), 3631–3636.
E = speye(m);
D = diff(E, d);
% C = chol(E + lambda * D' * D);
% z = C \ (C' \ y);

% Exact hat diagonal
H = inv(E + lambda * D' * D);
h = diag(H);
trace(H)
%% the third when m > 100
n = 100;
E1 = speye(n);
D1 = diff(E1, d);
lambda1 = lambda * (n / m) ^ (2 * d);
H1 = inv(E1 + lambda1 * D1' * D1);
% h1 = diag(H1)
trace(H1)

% In the command window: 
%   The first is  147.5582
%   The second is 73.8645
%   The third is  83.4209
%   -------------------------------------------------------------
%   And, in their paper, the free parameters of Whittaker smoother 
%   is 18.84 when λ=2, and it is 9.37 when λ=15.
