function freeEnergy = freeEnergy(phi,alphaL,alphaR,theta)
    syms E_syms
    DeltaL = cos(2*(theta-alphaL));
    DeltaR = cos(2*(theta-alphaR));
    etaL = acos(E_syms/abs(DeltaL));
    etaR = acos(E_syms/abs(DeltaR));
    if DeltaL < 0
        etaL = etaL + pi;
    elseif DeltaR < 0
        etaR = etaR + pi;
    end
    assume(E_syms,'real')
    E = vpasolve(cos(etaL + etaR)-cos(phi)==0,E_syms,[0 1],'random',true);
    freeEnergy = -2*log(2*cosh(E/2));
end