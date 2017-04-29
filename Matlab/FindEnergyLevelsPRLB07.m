disp(1)
syms k E xi real
nu = 40;
delta = 2;
Z = 0;
sgn = 1;

N = 2;
k_array = linspace(-1.5,1.5,N);
E1_array = zeros(size(k_array));
E2_array = zeros(size(k_array));

eta = 2/nu*sqrt(delta^2-E^2);
q_e = sqrt(1-k^2-1i*eta);
q_h = sqrt(1-k^2+1i*eta);
xi_e = sqrt(2*nu)*sgn*k;
xi_h = -xi_e;
gamma_e = delta./(E+1i*sqrt(delta^2-E^2));
gamma_h = delta./(E-1i*sqrt(delta^2-E^2));
a_e = -(nu/2+E);
a_h = -(nu/2-E);
U_e = xi^(-0.5)*kummerU(-a_e/2,-0.25,xi^2/2);
U_h = xi^(-0.5)*kummerU(-a_h/2,-0.25,xi^2/2);
dU_e = sqrt(2/nu)*diff(U_e,xi);
dU_h = sqrt(2/nu)*diff(U_h,xi);
U_e = subs(U_e,xi,xi_e);
U_h = subs(U_h,xi,xi_h);
dU_e = subs(dU_e,xi,xi_e);
dU_h = subs(dU_h,xi,xi_h);

matrix = [U_e 0 -gamma_e -gamma_h; 0 U_h -1 -1; dU_e-Z 0 -1i*q_e*gamma_e 1i*q_h*gamma_h; 0 dU_h-Z -1i*q_e 1i*q_h];
determinant = det(matrix);

disp(2)

for i=1:N
    tempDet = subs(determinant,k,k_array(i));
    re = real(tempDet);
    im = imag(tempDet);
disp(3)
    funReal = matlabFunction(re);
    funImag = matlabFunction(im);
disp(4)
    E1_array(i) = fzero(funReal,1.4);
disp(5)
    %E2_array(i) = fzero(funImag,1.1);
end

plot(k_array,E1_array);

