function [omega1hat, omega1] = genJointDist(polp, pollamb, params, options)

    arguments
        polp double
        pollamb double
        params struct
        options.maxiter double = 1000
        options.tol double = 1e-6
        options.printinterval double = 100
        options.printinfo logical = true
    end

    aP = params.aP;

    omega0hat = ones(params.np, params.na);
    omega0hat = omega0hat ./ (params.np*params.na);
    error = 10;
    iter = 0;
    
    while (error > options.tol) && (iter < options.maxiter)

        % update p dist
        omega1 = zeros(params.np, params.na);

        for pidx = 1:params.np
            for aidx = 1:params.na
                
                p1idx = polp(pidx, aidx);

                % non adjusters
                omega1(pidx, aidx) = omega1(pidx, aidx) + (1 - pollamb(pidx, aidx)) * omega0hat(pidx, aidx);
                
                % adjusters
                omega1(p1idx, aidx) = omega1(p1idx, aidx) + pollamb(pidx, aidx) * omega0hat(pidx, aidx);
            end
        end
        
        % update shock dist
        omega1hat = zeros(params.np, params.na);

        for pidx = 1:params.np
            for aidx = 1:params.na

                for a0idx = 1:params.na
                    omega1hat(pidx, aidx) = omega1hat(pidx, aidx) + omega1(pidx, a0idx) * aP(a0idx, aidx);
                end
           end
        end

        error = max(abs(log(1+omega1hat) - log(1+omega0hat)), [], "all");
        iter = iter + 1;
        omega0hat = omega1hat;

        if ((iter==1) || (mod(iter, options.printinterval) == 0)) & options.printinfo
            fprintf("Iterations: %d, Error %5.9f\n", iter, error)
        end

    end

    if options.printinfo
        fprintf("Final Iterations: %d, Final error: %5.9f\n", iter, error)
    end

end
