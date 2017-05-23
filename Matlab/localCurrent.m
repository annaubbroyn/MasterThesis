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
        
        %alphaR = alphaL;
        if alphaL == alphaR 
            alpha = alphaL;
            E = sqrt(cos(2*theta-2*alpha)^2+swave^2)*cos(phi/2-gamma/2);
            dE = sqrt(cos(2*theta-2*alpha)^2+swave^2)*sin(phi/2-gamma/2);
            I = 1/(1+swave^2)*tanh(E/2)*dE;
        elseif alphaL == pi/4 && alphaR == -pi/4
            E = sqrt(sin(2*theta)^2+swave^2)*sin(phi/2-gamma/2+atan(swave/sin(2*theta)));
            dE = sqrt(sin(2*theta)^2+swave^2)*cos(phi/2-gamma/2+atan(swave/sin(2*theta)));
            I = 1/(1+swave^2)*tanh(E/2)*dE;
        end
            %theta = 0;
            %I = cos(2*(theta-alpha))*sin(phi/2-gamma/2)*tanh(cos(2*(theta-alpha))*cos(phi/2-gamma/2)/2);
    end
end