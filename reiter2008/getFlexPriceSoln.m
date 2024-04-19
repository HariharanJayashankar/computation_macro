function [Yf, PDf, wf, Lf] = getFlexPriceSoln(p)

    Yf = 1; % this is normmalized


    % equation 16 in pset
    deltaf = (exp(p.agrid)).^(p.epsilon - 1) * p.aPstar;
    PDf = deltaf^(1/(p.epsilon - 1));

    % equation 14
    Lf = PDf * Yf;

    wf = ((p.epsilon - 1)/p.epsilon) * (1/PDf);



end
