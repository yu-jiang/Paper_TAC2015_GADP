function dx = mySusp(t,x,u,noise_on)
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m

% State matrices
A = [ 0 1 0 0;
    [-ks -bs ks bs]/mb ; ...
    0 0 0 1;
    [ks bs -ks-kt -bs]/mw];
B = [ 0; 10000/mb ; 0 ; -10000/mw];
B1 = [ 0; 0 ; 0 ; kt/mw];
x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);
if noise_on == 1
    e=0.005*sum(sin([38.1558   76.5517   79.5200   18.6873   48.9764   44.5586   64.6313   70.9365   75.4687   27.6025]*t));
    %e2=0.005*sum(sin([17.9703   15.5098  -33.7388  -38.1002   -0.1636   45.9744  -15.9614    8.5268  -27.6188   25.1267]*t));
else
    e =0;
    %e2 =0;
end
%dx = A*x+B*eval(u)/(1+x'*x)+ B*e;
if abs(t-120)<=0.01
    r= 0;
else
    r= 0;
end
dx = A*x+B*eval(u)+ B*0*e+B1*r;
end