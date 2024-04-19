function lout = labour_static(c, w, p)
    % uses static hh optimality to give
    % labour supply igven a consumption nad wages

    lout = (w/(c*p.zeta))^p.nu;

    

end
