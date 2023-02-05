%% Basic parameters of the setup

% All parameters and variables are named approximately as in my master's
% thesis. The parameter values chosen below reproduce the one of the graphs
% shown in Figure 4.4 of the thesis.

gamma = 1/pi;

r1 = 0.5; % inner radius of the nonscattering region
r2 = 0.8; % outer radius of the nonscattering region
b = 0.9; % 'anisotropy' coefficient (denoted by B in the thesis)
mua = 0.5; % absorption coefficient in the strongly scattering region
mus = 500; % scattering coefficient in the strongly scattering region
K = 1/(2*(mua + (1-b)*mus)); % diffusion coefficient in the strongly scattering region
muan = 0.25; % absorption coefficient in the nonscattering region

%% Parameters deduced from those above

delta = sqrt(mua/K); % parameter appearing in the arguments of the Bessel functions

% The two integrands needed for the evaluating the integrals defining g1, g2
% and g3. There are two small errors in the thesis that are fixed below: (i) in the integrals
% defining g1 and g2, the normalization by |S^0| is missing in the thesis, which
% correponds to division by 2. This corresponds to the extra multiplication by 0.5
% in the definition of fun12 below. (ii) r2 should appear in the
% exponential term of the integral defining g3, which is now accounted for
% in fun3.

fun12 = @(theta) ( 0.5*(r1^2*cos(theta) - r1*r2*(1+cos(theta).^2) + r2^2*cos(theta))./(r1^2 - 2*r1*r2*cos(theta) + r2^2).^(3/2) ).*exp(-muan*(r1^2 - 2*r1*r2*cos(theta) + r2^2).^(1/2));
fun3 = @(theta) (1 - cos(theta)).^(1/2).*exp(-muan*r2*(2*(1-cos(theta))).^(1/2));

% the corresponding upper limits for integration
b12 = acos(r1/r2);
b3 = 2*b12;

% evaluating the values of g1, g2 and g3 numerically
int12 = integral(fun12, 0, b12, 'AbsTol', 1e-12);
g1 = -2*r2*int12;
g2 = -2*r1*int12;

int3 = integral(fun3, 0, b3, 'AbsTol', 1e-12);
g3 = -1/(2*sqrt(2))*int3;

% the coefficients a1, a2, a3 and a4 according to (4.16)
a1 = (1 + g1*g2 - g3)/(1 - g1*g2 - g3);
a2 = 2*g1/(1 - g1*g2 - g3);
a3 = 2*g2/(1 - g1*g2 - g3);
a4 = (1 + g1*g2 + g3)/(1 - g1*g2 - g3);


%% Deducing the coefficients for the linear equations (4.17), (4.18) and (4.19)

A = zeros(3); % the system matrix

% The first variable is c2I, the second c2K and the third c1I

% The first row is (4.17). 
A(1,1) = 2*gamma*besseli(0,delta) + K*delta*besseli(1,delta);
A(1,2) = 2*gamma*besselk(0,delta) - K*delta*besselk(1,delta);
A(1,3) = 0;

% The second row is (4.18)
A(2,1) = 2*gamma*a2*besseli(0,delta*r2);
A(2,2) = 2*gamma*a2*besselk(0,delta*r2);
A(2,3) = 2*gamma*a1*besseli(0,delta*r1) + K*delta*besseli(1,delta*r1);

% The third row is (4.19). As you noted, r1 should be r2 at one point in
% (4.19):
A(3,1) = 2*gamma*a4*besseli(0,delta*r2) - K*delta*besseli(1,delta*r2);
A(3,2) = 2*gamma*a4*besselk(0,delta*r2) + K*delta*besselk(1,delta*r2);
A(3,3) = 2*gamma*a3*besseli(0,delta*r1);

% The vector on the right-hand side
y = [2; 0; 0];

% solving for the coefficients of the modified Bessel functions
c = A\y;

%% Plotting the energy fluency

% the area outside the nonscattering regions
rho1 = linspace(r2,1,100);
phi1 = c(1)*besseli(0,delta*rho1) + c(2)*besselk(0,delta*rho1);

% the area inside the nonscattering region
rho2 = linspace(0,r1,100);
phi2 = c(3)*besseli(0,delta*rho2);

figure(1)
plot(rho1,phi1);
hold on
plot(rho2,phi2);
hold off



