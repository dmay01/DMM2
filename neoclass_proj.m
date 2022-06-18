% solves the stochastic neoclassical growth model using projection methods

%% !!! install the CompEcon toolbox first (see slides)

%% define parameter values and functional forms
% initialize P
P = struct();
% define parameters
P = pars(P);
% define functions
P = funforms(P);

%% approximate the value function (B7)
% settings
P.nK = 4;      % number of polynomials used in K dimension (must be >=2)
P.nA = 5;      % number of polynomials used in A dimension (must be >=2)
P.nquad = 5;   % number of nodes for Gauss-Hermite quadrature
% setup the intervals and basis functions
K = [0.8,1.4];
A = [0.5,1.2];
fspace = fundefn('cheb',[P.nK,P.nA],[K(1),A(1)],[K(2),A(2)]);
% compute the Chebyshev nodes on the specified intervals
x = funnode(fspace);
% compute the projection
cstar = compute_proj(P,fspace,x);

%% check accuracy of the approximation (B8)
% 100 evenly spaced grid points in each dimension
[~,xcheck] = nodeunif([100,100],fspace.a,fspace.b);

% evaluate Euler residuals
P.nquad = 10;   % increase nodes for quadrature
rescheck = compute_residual(P,fspace,xcheck,cstar);

% maximum residual
maxerr = max(rescheck);
fprintf('max. residual:\t%6.2e \n', maxerr)

% minimum residual
minerr = min(rescheck);
fprintf('min. residual:\t%6.2e \n', minerr)

%% plot approximated value function B9
A = 0.5:0.5/100:1;
K = 1.110432704612086;
Vq = ones(size(A));
for i = 1:length(A)
    Vq(i) = initialguess_quad([K, A(i)]);
end
V = funeval(cstar,fspace, [K*ones(length(A),1),nodeunif(length(A),0.5,1)]);
plot(A,Vq);
hold on
plot(A,(V));
xlabel('productivity A');
ylabel('value V');
title('approximated value function Vhat(K,A) vs. initial quadratic guess Vq(K,A)');
legend('Vq','Vhat')

%% welfare impact of approximation B10
sh = ones(size(A));
for i = 1:length(A)
    sh(i) = fsolve(@(s) P.u(1-s)/(1-P.beta) + (1-s)*P.du(1-s)*Vq(i) -V(i),0.5);
end
plot(A,sh);
xlabel('productivity A');
ylabel('share s');
title('impact of quadratic approximation as share of consumption variation');