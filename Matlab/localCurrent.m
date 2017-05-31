function I = localCurrent(x,y,theta,l,phi,option,lambda,system,alphaL,alphaR,h,swave)
    L = 100;
    gamma = Gamma(x,y,theta,l,option,lambda);
    if strcmp(system,'SNS')
        I = sin(phi/2-gamma/2)*tanh(cos(phi/2-gamma/2)/2);
        xi =  h*L*cos(theta);
    elseif strcmp(system,'SFS')
        %I = sin(phi/2-gamma/2-sigma*xi/2)*tanh(cos(phi/2-gamma/2-sigma*xi/2)/2);  
        Ip = sin(phi/2-gamma/2+xi/2)*tanh(0.5*cos(phi/2-gamma/2+xi/2));
        Im = sin(phi/2-gamma/2-xi/2)*tanh(0.5*cos(phi/2-gamma/2-xi/2));
        I = 0.5*(Ip+Im);
        %I = 0.5*sin(phi/2-gamma/2+xi/2)/(cosh(cos(phi/2-gamma/2+xi/2))*cosh(cos(phi/2-gamma/2-xi/2)));
    elseif strcmp(system,'SO')
        
    elseif strcmp(system,'Dwave')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        deltaMax = sqrt(1+swave^2);
        d = swave;
        b = 2*theta;
        cR = 2*alphaR;
        cL = 2*alphaL;
        chiRp = atan(d/cos(b-cR)) + pi*(cos(b-cR)<0);
        chiRm = atan(d/cos(b+cR)) + pi*(cos(b+cR)<0);
        chiLp = atan(d/cos(b-cL)) + pi*(cos(b-cL)<0);
        chiLm = atan(d/cos(b+cL)) + pi*(cos(b+cL)<0);
        A_p = sqrt(cos(b+cL)^2+d^2)*sqrt(cos(b-cR)^2+d^2);
        A_m = sqrt(cos(b-cL)^2+d^2)*sqrt(cos(b+cR)^2+d^2);
        B_p = cos(b+cL)^2+cos(b-cR)^2+2*d^2;
        B_m = cos(b-cL)^2+cos(b+cR)^2+2*d^2;
        denom_p = sqrt(B_p-2*A_p*cos(chiRp+chiLm+phi-gamma));
        denom_m = sqrt(B_m-2*A_m*cos(chiRm+chiLp-phi+gamma));
        numer_p = A_p*sin(chiRp+chiLm+phi-gamma);
        numer_m = A_m*sin(chiRm+chiLp-phi+gamma);
        if denom_p > 0
            ap = (1/deltaMax)*numer_p/denom_p;
            dap = (1/deltaMax)*(A_p^2*(cos(2*(chiRp+chiLm+phi-gamma))+3)-2*A_p*B_p*cos(chiRp+chiLm+phi-gamma))/(2*(B_p-2*A_p*cos(chiRp+chiLm+phi-gamma))^(3/2));
        elseif abs(numer_p)<10^(-10)
            ap = 0;
            dap = 0;
        else
            disp('Dividing on zero in ap')
            disp('denom_p')
            disp(denom_p)
            disp('numer_p')
            disp(numer_p)
        end
        if denom_m > 0
            am = (1/deltaMax)*numer_m/denom_m;
            dam = -(1/deltaMax)*(A_m^2*(cos(2*(chiRm+chiLp-phi+gamma))+3)-2*A_m*B_m*cos(chiRm+chiLp-phi+gamma))/(2*(B_m-2*A_m*cos(chiRm+chiLp-phi+gamma))^(3/2));
        elseif abs(numer_m)<10^(-10)
            am = 0;
            dam = 0;
        else
            disp('Dividing on zero in am')
            disp('denom_m')
            disp(denom_m)
            disp('numer_m')
            disp(numer_m)
        end
        
        
        Ip = -tanh(ap/2)*dap;
        Im = -tanh(am/2)*dam;
        I = Ip + Im;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           E = abs(sin(2*theta))*cos(phi/2-gamma/2);
%           dE = -0.5*abs(sin(2*theta))*sin(phi/2-gamma/2);
%           I = -2*tanh(E/2)*dE;
%         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%         %h = 0.001;
%         %gamma = 0;
%         %alphaL = alphaR;
%         %f1 = freeEnergy(phi+h-gamma,alphaL,alphaR,theta);
%         %f2 = freeEnergy(phi-gamma,alphaL,alphaR,theta);
%         %I = (f1-f2)/h;
%         
%         %alphaR = alphaL;
%         if alphaL == alphaR 
%             alpha = alphaL;
%             E = sqrt(cos(2*theta-2*alpha)^2+swave^2)*cos(phi/2-gamma/2);
%             dE = sqrt(cos(2*theta-2*alpha)^2+swave^2)*sin(phi/2-gamma/2);
%             I = 1/(1+swave^2)*tanh(E/2)*dE;
%         elseif alphaL == pi/4 && alphaR == -pi/4
%             E = sqrt(sin(2*theta)^2+swave^2)*sin(phi/2-gamma/2+atan(swave/sin(2*theta)));
%             dE = sqrt(sin(2*theta)^2+swave^2)*cos(phi/2-gamma/2+atan(swave/sin(2*theta)));
%             I = 1/(1+swave^2)*tanh(E/2)*dE;
%         end
%             %theta = 0;
%             %I = cos(2*(theta-alpha))*sin(phi/2-gamma/2)*tanh(cos(2*(theta-alpha))*cos(phi/2-gamma/2)/2);
    end
end