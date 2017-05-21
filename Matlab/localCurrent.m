function I = localCurrent(x,y,theta,l,phi,option,lambda,system,alphaL,alphaR,h,swave)
    sigma = 1;
    L = 100;
    factor = h*L;
    xi =  factor*(1+tan(theta)^2)*cos(theta);
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
        %h = 0.001;
        %gamma = 0;
        %alphaL = alphaR;
        %f1 = freeEnergy(phi+h-gamma,alphaL,alphaR,theta);
        %f2 = freeEnergy(phi-gamma,alphaL,alphaR,theta);
        %I = (f1-f2)/h;
        alpha = alphaL;
        if alphaL==alphaR
            E = sqrt(cos(2*(theta-alpha))^2+swave^2)*cos(phi/2-gamma/2);
            dE = -0.5*sqrt(cos(2*(theta-alpha))^2+swave^2)*sin(phi/2-gamma/2);
        elseif abs(alphaL-alphaR)==pi/2
            chi = atan(swave/cos(2*(theta-alpha)));
            E = sqrt(cos(2*(theta-alpha))^2+swave^2)*cos(phi/2-gamma/2+chi);
            dE = -0.5*sqrt(cos(2*(theta-alpha))^2+swave^2)*sin(phi/2-gamma/2+chi);
        elseif abs(alpha-alphaR)==pi/4
            f=sqrt(cos(2*(theta-alpha))^2*sin(2*(theta-alpha))+swave^2+swave^4);
            E = f*sin(phi-gamma+chiL-chiR)/sqrt(2*swave^2+1-f*cos(phi-gamma+chiL-chiR));
            dE = f*(cos(phi-gamma+chiL-chiR)/sqrt(2*swave^2+1-sgn*f*cos(phi-gamma+chiL-chiR)
            %theta = 0;
            %I = cos(2*(theta-alpha))*sin(phi/2-gamma/2)*tanh(cos(2*(theta-alpha))*cos(phi/2-gamma/2)/2);
    end
end