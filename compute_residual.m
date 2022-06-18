function res = compute_residual(P,fspace,x,coef)
    % compute the Euler residuals at the gridpoints x when coef is used as 
    % the coefficient vector in our polynomial basis fspace
    
    % obtain quadrature nodes and weights
    [eps,w]=qnwnorm(P.nquad,0,P.sigma2);  

    % load the grid points
    K = x{1};
    A = x{2};
    % get the dimensions
    mK = numel(x{1});
    mA = numel(x{2});
    
    % initialize residuals
    R = nan(mK,mA);

    % loop through grid points in K dimension
    for i=1:mK
        % loop through grid points in A dimension
        for j=1:mA
            % evaluate V(K_i,A_j)
            V = funeval(coef,fspace,[K(i),A(j)]);
            % evaluate VK(K_i,A_j)
            VK = funeval(coef,fspace,[K(i),A(j)], [1,0]);
            % express C(K,A)
            C = ((1-P.delta+A(j)*P.FK(K(i),P.l))/VK)^(1/P.eta);
            % evaluate K'(K_i,A_j)
            Knext = (1-P.delta)*K(i)+A(j)*P.F(K(i),P.l)-C;
            % evaluate A'(K_i,A_j) at the quadrature nodes (eps)
            Anext = P.g(A(j),eps);     
            % evaluate V'(K_i,A_j)=V(K'(K_i,A_j),A'(K_i,A_j)) at the quadrature nodes
            Vnext = funeval(coef,fspace,[Knext*ones(size(Anext)),Anext]); 
            % approximate the integral
            I = w'*Vnext;
            % evaluate the Bellman residual
            R(i,j) = P.u(C) + P.beta*I - V;
        end
    end
    % turn into a vector
    res = R(:);

end
