F=[9/8 -1 -1/2; 0 0 0];
G=[0;1];
R=1;
Q=diag([1,1,0]);
%K=[6-3/4 -3 -1];
K=[6 -3 0];
%p1=5/2;
%p2=-1;%-1;-1 1/2];

%a=0
%b=1/4
cvx_begin sdp
cvx_precision high
variable P(2,2) symmetric
variable a
variable b
P>=0;
a>=0;
W=[2*P(1,1) P(1,2)+P(2,1) a;P(1,2)+P(2,1) 2*P(2,2) 0];
H=-(1/2*W'*(F+G*K)+1/2*(F+G*K)'*W+Q+K'*R*K);
HH=H-[0 0 b; 0 0 0; b 0 0];
b>=0;
HH>=0;
cvx_end
% %%
% syms x1 x2 sigma real
% V=1/2*x1*x1+1/2*(x2-2*x1)^2
% Vx=[diff(V,x1)  diff(V,x2)];
% Fc=F+G*K;
% Fcx=Fc*[x1 x2 sigma]'
% LfV=simple(Vx*Fcx)
% 
% %%
% 
% P=[5/2 -1;-1 1/2];
% W=[2*P(1,1) P(1,2)+P(2,1) a;P(1,2)+P(2,1) 2*P(2,2) 0];
% H=-(1/2*W'*(F+G*K)+1/2*(F+G*K)'*W);%+Q+K'*R*K)
% HH=H-[0 0 b; 0 0 0; b 0 0]
