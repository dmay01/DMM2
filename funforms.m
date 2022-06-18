function P = funforms(P)
% defines production and utility function etc.
% inputs: parameters P
% output: P updated with several functions

% production function
P.F  = @(K,L) K.^P.alpha.*L.^(1-P.alpha);
% marginal products
P.FK = @(K,L) P.alpha*K.^(P.alpha-1).*L.^(1-P.alpha);
P.FL = @(K,L) (1-P.alpha)*K.^P.alpha.*L.^(-P.alpha);

% utility function
if P.eta == 1
    P.u = @(c) log(c);
else
    P.u = @(c) (c.^(1-P.eta)-1)/(1-P.eta);
end

% derivate of utility
P.du = @(c) c.^(-P.eta);

% producitivity process
P.g = @(A,eps) P.rho*A + (1-P.rho)*P.Abar +P.sigma* eps;

end