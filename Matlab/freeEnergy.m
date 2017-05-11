function freeEnergy = freeEnergy(phi,alphaL,alphaR,theta)
    syms E_syms
    etaL = acos(E_syms/cos(2*(theta-alphaL)));
    etaR = acos(E_syms/cos(2*(theta-alphaR)));
    assume(E_syms,'real')
    E = vpasolve(cos(etaL + etaR)-cos(phi)==0,E_syms,[0 1],'random',true);
    freeEnergy = -2*log(2*cosh(E/2));
end