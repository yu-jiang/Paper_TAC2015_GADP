function dx = mySusp0(t,x)
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

B = [ 0; 0 ; 0 ; kt/mw];

if abs(t-120)<=0.01
    r= 0;
else
    r= 0;
end
%
dx = A*x + B*r;

end