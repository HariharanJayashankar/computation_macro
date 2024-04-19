function pout = flexsoln(w, params)

    % get flexi price solution (which is the static solution)
    % w = wage
    % outputs a price vector of same sizee as the shock grid

    pout = (params.epsilon/(params.epsilon - 1)) .* w ./ (exp(params.agrid));
    
end
