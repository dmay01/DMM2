// First we define the variables (exogenous and endogenous) and the parameters 


var C K A V;
varexo e;

parameters beta eta alpha delta rho A_line sigma l;
beta    = 0.9;
eta     = 2;
alpha   = 0.25;
delta   = 0.12;
rho     = 0.9; 
A_line  = 1;
sigma   = 0;
l = 1;



//Declaring the model
model;
    A       = rho*A(-1) + (1-rho)*A_line + sigma*e;
    K        = (1-delta)*K(-1) + A*K(-1)^alpha*l^(1-alpha) - C;
    C^(-eta) = beta*(1-delta+alpha*A(+1)*K^(alpha-1)*l^(1-alpha))*C(+1)^(-eta);
    V           = ((C^(1-eta)-1)/(1-eta)) + beta*V(+1);
 
end;
// Steady 

initval;
    C = 1;
    K = 3;
    A = 1;
    V = 1;
end;

steady;

check;

// specify shocks

shocks;
    var e; stderr 1;
end;

// simulate

// B3
stoch_simul(order=1,nograph,nofunctions,nocorr,ar=0);
set_param_value('sigma',0.02);
stoch_simul(order=1);

// B4
set_param_value('sigma',0);
stoch_simul(order=2);
set_param_value('sigma',0.02);
stoch_simul(order=2);

// B5
set_param_value('eta',5);
set_param_value('sigma',0);
stoch_simul(order=2);
set_param_value('sigma',0.02);
stoch_simul(order=2);

// B6
set_param_value('eta',2);
set_param_value('rho',0.95);
set_param_value('sigma',0);
stoch_simul(order=2);
set_param_value('sigma',0.02);
stoch_simul(order=2);