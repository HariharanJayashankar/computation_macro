function [vout, polpout, pollambdaout, error, iter] = viterFirm(agg, params, options)
    % value function iterations
    % v(p, a)
    % where p is firms last price
    % a is shock value
    % Gamma is aggregate state
    % assume the steady state equilibrium we have is s.t.
    % aggregate dist doesnt change, therefore all the firm
    % needs is aggregate Y and wages, and not the whole joint distribution
    % of p and a
    % agg is a struct with fields w and Y


    arguments
        agg struct
        params struct
        options.maxiter double = 10000
        options.tol double = 1e-6
        options.printinterval double = 1000
        options.printinfo logical = true
    end
    

    % create profit matrix for future use
    pi_mat = zeros(params.np, params.na);
    for pidx=1:params.np
        for aidx=1:params.na
            pval = params.pgrid(pidx);
            aval = params.agrid(aidx);
            pi_mat(pidx, aidx) = (pval^(1-params.epsilon) - pval^(-params.epsilon)*(agg.w/exp(aval))) * agg.Y;
        end
    end

    % initial values for v
    vadjust = zeros(params.np, params.na);
    vnoadjust = zeros(params.np, params.na);

    error = 10;
    iter = 0;

    while (error > options.tol) && (iter < options.maxiter)

        v0 = max(vadjust, vnoadjust);
        expectedV = v0 * params.aP';

        vnoadjust = pi_mat + params.beta * expectedV;
        valadjust_tmp = pi_mat - params.kappa + params.beta * expectedV; 

        % valadjust_tmp is np * na - need to get max for any given a
        [vadjust, polpout_tmp] = max(valadjust_tmp, [], 1);
        vadjust = repmat(vadjust, params.np, 1);

        v1 = max(vadjust, vnoadjust);
        error = max(abs(v1 - v0), [], 'all');
        iter = iter + 1;

        if (iter == 1 || mod(iter, options.printinterval) == 0) & options.printinfo
            fprintf("Iterations: %d, Error: %5.9f\n", iter, error)
        end
        

    end

    if iter == options.maxiter
        fprintf("Warning: max iterations reached. Problem may have not converged\n")
    end

    if options.printinfo
        fprintf("Firm Viter: Final iterations %d, Final error %5.9f\n", iter, error)
    end

    polpout = repmat(polpout_tmp, params.np, 1); 
    pollambdaout = vadjust > vnoadjust;
    vout = v1;

end
