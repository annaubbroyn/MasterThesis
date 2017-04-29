X;
k;
beta=0;
dhw=2;
syms a b z xi E
g = (sqrt(1+(dhw/k^2))+1)/(sqrt(1+(dhw/k^2))-1);
l = sqrt(abs(X)/k);
nu = 40;
alpha = -l^2/2*(nu/2+E);
W = whittakerW(a,b,z);
dW = diff(W,z);
W = subs(W,a,-alpha/2);
W = subs(W,b,-1/4);
W = subs(W,z,xi^2/2);
dW = subs(dW,a,-alpha/2);
dW = subs(dW,b,-1/4);
dW = subs(dW,z,xi^2/2);
U = 1/sqrt(xi)*2^(-alpha/2)*W;
dU = -0.5*xi^(-3/2)*2^(-alpha/2)*W+1/sqrt(xi)*2^(-alpha/2)*dW;
UzeroP = subs(U,xi,sqrt(2)*X/l);
dUzeroP = subs(dU,xi,sqrt(2)*X/l);
UzeroM = subs(UzeroM,xi,-sqrt(2)*X/l);
dUzeroM = subs(dUzeroM,xi,-sqrt(2)*X/l);
UzeroM = subs(U,E,-E);
dUzeroM = subs(dU,E,-E);

matrix = [UzeroP 0 -1 g; 0 UzeroM -g -1; dUzeroP 0 -1i*k+beta (1i*k+beta)*g; 0 dUzeroM (-1i*k+beta)*g 1i*k+beta];

det(matrix);

