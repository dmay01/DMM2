/*
 * Stochastic neoclassical growth model with endogenous labour supply.
 * author: Daniel May
 */

// I. list variabes and paramters

var C K A L Y I;
varexo e;

parameters beta alpha delta rho theta psi;
beta = 0.995;
alpha = 0.33;
delta = 0.012;
rho   = 0.95;
theta = 0.873074970384159;
psi= 0.5;

// II. specify the model

model;  
    1/C = beta*(1-delta+alpha*A(+1)*K^(alpha-1)*L(+1)^(1-alpha))/C(+1);
    K        = (1-delta)*K(-1) + A* K(-1)^alpha*L^(1-alpha) - C;
    log(A)   = rho*log(A(-1)) + e;
    L^psi = (A*(1-alpha)*K(-1)^alpha*L^(-alpha))/(theta * C);
    Y  = A*K(-1)^alpha*L^(1-alpha);
    I = Y - C;
end;

// III. solve for the steady state

initval;
C = 1;
K = 3;
A = 1;
L = 1;
I = 1;
Y = 1;
end;

// show steady state
steady;

// check rank condition
check;

// IV. specify shocks

shocks;
    var e; stderr 0.009;
end;

// V. simulate and plot impulse responses

stoch_simul(order=1,IRF=200, hp_filter=1600);
