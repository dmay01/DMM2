function cstar = compute_proj(P,fspace,x)
    % compute the projection on the polynomial basis defined in fspace,
    % using the gridpoints specified by x

    % fit the guessed function (in initialguess.m) with Chebyshev interpolation
    % to get a representation in terms of our Chebyshev polynomial basis
    coef0 = funfitf(fspace,"initialguess_quad");
    % this is the initial guess for the coefficient vector  

       
    % do the projection
    cstar = fsolve(@(c) compute_residual(P,fspace,x,c), coef0);
    
end