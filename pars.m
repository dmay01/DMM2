function P = pars(P)
% defines the model parameters
% inputs: structure P
% output: P updated with parameter values

% discount factor
P.beta = 0.9;
% relative risk aversion
P.eta = 2;
% labor supply (normalization)
P.l = 1;
% output elasticity of capital
P.alpha = 0.25;
% depreciation rate
P.delta = 0.12;
% persistence of productivity process
P.rho = 0.9;
% constant for productivity shocks
P.sigma = 0.02;
% var of productivity shocks
P.sigma2 = 1;
% steady state productivity
P.Abar = 1;

end