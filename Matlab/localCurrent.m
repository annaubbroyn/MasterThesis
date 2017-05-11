function I = localCurrent(x,y,theta,l,phi,option,lambda,system,alphaL,alphaR)
    sigma = 1;
    L = 100;
    xi =  cos(theta)/(2*l^2*L);
    gamma = Gamma(x,y,theta,l,option,lambda);
    if strcmp(system,'SNS')
        I = sin(phi/2-gamma/2)*tanh(cos(phi/2-gamma/2)/2);
    elseif strcmp(system,'SFS')
        %I = sin(phi/2-gamma/2-sigma*xi/2)*tanh(cos(phi/2-gamma/2-sigma*xi/2)/2);  
        Ip = sin(phi/2-gamma/2+xi/2)*tanh(0.5*cos(phi/2-gamma/2+xi/2));
        Im = sin(phi/2-gamma/2-xi/2)*tanh(0.5*cos(phi/2-gamma/2-xi/2));
        I = 0.5*(Ip+Im);
        %I = 0.5*sin(phi/2-gamma/2+xi/2)/(cosh(cos(phi/2-gamma/2+xi/2))*cosh(cos(phi/2-gamma/2-xi/2)));
    elseif strcmp(system,'Dwave')
        h = 0.001;
        gamma = 0;
        f1 = freeEnergy(phi+h-gamma,alphaL,alphaR,theta);
        f2 = freeEnergy(phi-gamma,alphaL,alphaR,theta);
        I = (f1-f2)/h;
        %alpha = alphaL;
        %I = cos(2*(theta-alpha))*sin(phi/2-gamma/2)*tanh(cos(2*(theta-alpha))*cos(phi/2-gamma/2)/2);
    end
end