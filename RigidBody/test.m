%function test()
clear all; % echo on;
global K0 I1 I2 I3 u


% Inertial

syms x1 x2 x3 real

I1 = 3;
I2 = 1;
I3 = 0.5;

F = [  0   0   0   0     0     0       0  (I2-I3)/I1 zeros(1, 11);
    0   0   0   0     0     0     (I3-I1)/I2      0 zeros(1, 11);
    0   0   0   0  (I1-I2)/I3     0       0  zeros(1, 12); ];


Q = diag([1 1 1 zeros(1, 19 - 3)]);


G = diag([I1 I2 I3]);

z =[ x1, x2, x3, x1^2, x1*x2, x2^2, x1*x3, x2*x3, x3^2, x1^3, x1^2*x2, ...
    x1*x2^2, x2^3, x1^2*x3, x1*x2*x3, x2^2*x3, x1*x3^2, x2*x3^2, x3^3]';

K0 = zeros(3,19);
K0(1,1) = -1;
K0(2,2) = -1;
K0(3,3) = -1;
K0(1,10) = -1;
K0(2,13) = -1;
K0(3,19) = -1;

K=K0;

f = F*z + G*K0*z;

vars = [x1; x2; x3];
prog = sosprogram(vars);

[prog,V] = sospolyvar(prog,  monomials([x1; x2; x3], 2:4), 'wscoeff');

prog = sosineq(prog,V);

expr = -[diff(V,x1) diff(V,x2) diff(V,x3)]*f-z'*Q*z;%% ;

prog = sosineq(prog,expr);

prog = sossolve(prog);

V0 = sosgetsol(prog,V)
Vold =V0;


uold = K0*z;

%
VEC =  monomials([x1; x2; x3], 2:4);
c = int(int(int(VEC,-1,1),-1,1),-1,1);



% x0 = [2;-1;3];
% [t, y] = ode45(@RBsys, [0,10], x0);
% plot(t,y(:,1),t,y(:,2),t,y(:,3) )



%for j = 1: 1

clear prog V

vars = [x1; x2; x3];
prog = sosprogram(vars);

[prog,V] = sospolyvar(prog, VEC, 'wscoeff');

prog = sosineq(prog,Vold - V);

expr = -[diff(V,x1) diff(V,x2) diff(V,x3)]*f-z'*Q*z-uold'*uold;%% ;

prog = sosineq(prog,expr);

prog = sossolve(prog);

myObj =[coeff_1 coeff_2 coeff_3 coeff_4 coeff_5 coeff_6 ...
    coeff_7 coeff_8 coeff_9 coeff_10 coeff_11 coeff_12 ...
    coeff_13 coeff_14 coeff_15 coeff_16 coeff_17 coeff_18 ...
    coeff_19 coeff_20 coeff_21 coeff_22 coeff_23 coeff_24 ...
    coeff_25 coeff_26 coeff_27 coeff_28 coeff_29 coeff_30 ...
    coeff_31]  *c;

prog = sossetobj(prog, myObj);

V = sosgetsol(prog,V)
Vold =V;




u = - 0.5*G'*[diff(V,x1) diff(V,x2) diff(V,x3)]' ;

f = F*z + G*u;
%end

