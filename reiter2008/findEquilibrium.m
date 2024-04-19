function [out,V, polp, pollamb, omegahat, omega, Y0, Pi, Ld] = findEquilibrium(w, p)

    % gettngg value funcion
    Y0 = 1;
    agg.w = w;
    errorY = 10;
    tolY = 1e-3;
    max_yiter = 100;
    iterY = 0;
    updaterate = 0.3;


    while iterY < max_yiter && errorY > tolY


        agg.Y = Y0;
        [V, polp, pollamb, ~, ~] = viterFirm(agg, p, 'maxiter', 10000, 'tol', 1e-6, 'printinfo', false);

        %% get joint distribution of prices and shocks
        [omegahat, omega] = genJointDist(polp, pollamb, p, 'printinfo', false);

        % get implied aggregate Y
        Yimplied = sum((repmat(p.pgrid', 1, p.na) .^ (-p.epsilon) .* Y0) .^ omega, 'all');

        errorY = abs(Yimplied - Y0);
        Y1 = updaterate * Y0 + (1-updaterate)*Yimplied;
        Y0 = Y1;

        iterY = iterY + 1;

    end

    if iterY == max_yiter
        fprintf("Warning, maximum iterations reached for output. May not have converged\n")
    end

    % get profits to give HH
    % get aggregate fixed cost payments
    % get labour demand
    % integrate
    Pi = sum((p.pgrid' .^ (1-p.epsilon) - p.pgrid' .^ (-p.epsilon) .* (w./exp(p.agrid))) .* omega, 'all' );
    Ld =  sum(p.pgrid' .^ (-p.epsilon) .* exp(-p.agrid) .* Y0 .* omega, 'all');

    l_consumer = @(x) w/(w*x + Pi) - p.zeta*x^(1/p.nu);
    Ls = fsolve(l_consumer, 1);


    labmkterror = Ld - Ls;

    % error = abs(labmkterror);
    error = labmkterror;
    out = error;

end
