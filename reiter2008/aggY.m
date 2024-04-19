function Y = aggY(omega, p, agg)

    % given individual policy funcitons

    yj = zeros(p.np, p.na);
    for pidx = 1:p.np
        for aidx = 1:p.na
            pval = p.pgrid(pidx);
            yj(pidx, aidx) = pval^(-p.epsilon)*agg.Y;
        end
    end

    Y = 0;
    for pidx = 1:p.np
        for aidx = 1:p.na
            Y = Y + yj(pidx, aidx)^((p.epsilon-1)/p.epsilon)*omega(pidx, aidx);
        end
    end
    Y = Y^(p.epsilon/(p.epsilon-1));
end

