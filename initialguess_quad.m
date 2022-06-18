% for the given parameters, the second order perturbation yields the
% following quadratic function

function fval = initialguess_quad(x)
    % x = [K,A] where K and A are vectors
    fval = -1.210042 + 1.392454*(x(:,1)-1.11043) + 0.135417/0.02*(x(:,2)-1) ...
                     - 0.548820*(x(:,1)-1.11043).^2 - 0.032524/0.02*(x(:,1)-1.11043).*(x(:,2)-1) ...
                     - 0.001707/0.02^2.*(x(:,2)-1).^2;
end