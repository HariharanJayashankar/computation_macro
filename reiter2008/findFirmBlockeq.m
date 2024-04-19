function [Yerror, V, polp, pollamb, omegahat, omega] = findFirmBlockeq(Y, w, p)

    % within firm block get equilibirum aggregate Y

    % gettngg value funcion
    agg.Y = Y;
    agg.w = w;
    [V, polp, pollamb, ~, ~] = viterFirm(agg, p, 'maxiter', 10000, 'tol', 1e-6, 'printinfo', 1000000);

    %% get joint distribution of prices and shocks
    [omegahat, omega] = genJointDist(polp, pollamb, p);

    Yimplied = aggY(omega, p, agg);

    Yerror = abs(Yimplied - agg.Y);

end
