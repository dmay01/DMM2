close all;
clearvars;

%% parameters
P = struct();

P.beta = 0.995;
P.psi= 0.5;
P.delta = 0.012;
P.alpha = 0.33;
P.rho = 0.95;
P.sigma2 = 0.009^2;


%% functions
% production function
P.F  = @(K,L) K.^P.alpha.*L.^(1-P.alpha);
% marginal products
P.FK = @(K,L) P.alpha*K.^(P.alpha-1).*L.^(1-P.alpha);
P.FL = @(K,L) (1-P.alpha)*K.^P.alpha.*L.^(-P.alpha);

%% A4. calibrate theta
P.Lss = 1;
P.Ass = fzero(@(A) log(A)-P.rho*log(A),[0.01,1]);
P.Kss = fzero(@(K) 1 - P.beta*(1-P.delta + P.Ass*P.FK(K,P.Lss)),[0.001,100]);
P.Css = P.Ass*P.F(P.Kss, P.Lss) - P.delta*P.Kss ;
P.theta = (P.Ass*P.FL(P.Kss, P.Lss))/(P.Css*P.Lss.^P.psi);
disp(P.Ass);
disp(P.Kss);
disp(P.Css);
disp(P.theta);