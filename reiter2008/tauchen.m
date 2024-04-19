function [P, xgrid] = tauchen(m, rho, sigma, N)

    % Approximating process x_t = rho * x_{t-1} + eps
    % where eps ~ N(0, sigma^2)
    % m : a scalar which defines how many sd's away from the mean
    %      we want to go
    % rho : ar1 process
    % mu : Mean of error process. Scalar
    % sigma : sd of error process. Scalar
    % N : number of grid points (default is 1000)

    arguments
        m double
        rho double
        sigma double
        N double = 1000
    end


    sigmax = (sigma^2/(1 - rho^2))^(1/2);
    xN = m*sigmax;
    x1 = -xN;
    xgrid = linspace(x1, xN, N);
    w = xgrid(2) - xgrid(1);

    P = zeros(N, N);

    for i=1:N
        for j=1:N

            xi = xgrid(i);
            xj = xgrid(j);

            if j == 1
                P(i, j) = normcdf((xj - rho*xi + w/2)/sigma );
            elseif j == N
                P(i, j) = 1 - normcdf((xj - rho*xi - w/2)/sigma);
            else
                P(i, j) = normcdf((xj - rho*xi + w/2)/sigma) - normcdf((xj - rho*xi - w/2)/sigma);
            end

        end
    end

end
