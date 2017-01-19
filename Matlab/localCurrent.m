function I = localCurrentCircle(x,y,theta,l,phi,option,lambda)
    gamma = Gamma(x,y,theta,l,option,lambda);
    I = sin(phi/2-gamma/2)*tanh(cos(phi/2-gamma/2)/2);
    %I = 0.25*sin(gamma);
%     I = 0;
%     for i = 1:2:9
%         I = I+(-1)^i*gamma^(2*i+1)/factorial(2*i+1);
%     end
end