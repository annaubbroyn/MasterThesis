function gamma = Gamma(x,y,theta,l,option,lambda)
    a = y-x*tan(theta);
    b = tan(theta)/2;
    coeff = 1/(l^2*b);
    factor = 1;
    if option == 1
        %Constant B-field: B = B\hat(z)
        gamma = 2/l^2*a;
        %gamma = factor*coeff*2*a*b;
    elseif option == 2
        phi = 0;
        gamma = -factor*0.5*lambda^2/(l^2*b*pi^2)*sin(2*pi*b/lambda)*cos(2*pi*a/lambda + phi);
    elseif option == 3
        phi = pi/2;
        gamma = -factor*0.5*lambda^2/(l^2*b*pi^2)*sin(2*pi*b/lambda)*cos(2*pi*a/lambda + phi);
    elseif option == 4
        phi = 0;
        gamma = factor*2*lambda/(pi*l^2)*(a*sin(pi/lambda)*sin(phi)+b*(lambda/pi*sin(pi/lambda)-cos(pi/lambda))*cos(phi));
   elseif option == 5
        phi = pi/2;
        gamma = factor*2*lambda/(pi*l^2)*(a*sin(pi/lambda)*sin(phi)+b*(lambda/pi*sin(pi/lambda)-cos(pi/lambda))*cos(phi));
    end
end
% 
%     
%     elseif option == 2
%         %B = B*sin(omega*y)\hat(z)
%         omega = 0.525;%0.25;%(2)*pi/b;%pi/b+pi/4;
%         %gamma = 4/(l^2*omega^2*tan(theta))*cos(omega*a)*sin(omega*b);
%         gamma = coeff*(-2/omega^2)*cos(omega*a)*sin(omega*b);
%     elseif option == 3
%         %B = B*y\hat(z)
%         gamma = (coeff/6)*((a+b)^3-(a-b)^3);
%     elseif option == 4
%         %as option 5 but with tan(theta)
%         omega = 0.5;
%         lambda = 2*pi/omega;
%         phi = 0;
%         gamma = -25*2*a/l^2*sin(omega*b)/(omega*b)*cos(omega*a)/(omega*a);
%         %gamma = -lambda/(l^2*pi)*sin(pi*tan(theta)/lambda)/(pi*tan(theta)/lambda)*cos(2*pi*y/lambda + phi-2*pi*x*tan(theta)/lambda);
%     elseif option == 5
%         %with tan(theta) = 0
%         lambda = 4*pi;
%         phi = 0;
%         gamma = -lambda/(l^2*pi)*cos(2*pi*y/lambda + phi);
%     elseif option == 6
%         %with y*cos(x...)
%         lambda = 1/0.5;
%         gamma = (2*a*lambda/(pi*l^2))*sin(pi/lambda);
%     end