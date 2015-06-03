clear all; % echo on;

global K0 I1 I2 I3
% Define the system parameters
I1 = 3;
I2 = 1;
I3 = 0.5;
F = [zeros(1,7)  (I2-I3)/I1 zeros(1, 11);
    zeros(1,6)  (I3-I1)/I2 zeros(1, 12);
    zeros(1,4)  (I1-I2)/I3 zeros(1, 14); ];
G = diag([I1 I2 I3]);
Q = diag([10 10 10 zeros(1, 19 - 3)]);
R = 1;

K0 = zeros(3,19);  % Initial stabilizing controller
K0(1,1) = -1;
K0(2,2) = -1;
K0(3,3) = -1;
K0(1,10) = -1;
K0(2,13) = -1;
K0(3,19) = -1;

K = K0;



%Ksave = zeros(3,19,10);
Ksave(:,:,1)=K0;
%% Obtain the initial lyapunov function
syms x1 x2 x3 real
z =[ x1, x2, x3, x1^2, x1*x2, x2^2, x1*x3, x2*x3, x3^2, x1^3, x1^2*x2, ...
    x1*x2^2, x2^3, x1^2*x3, x1*x2*x3, x2^2*x3, x1*x3^2, x2*x3^2, x3^3]';
vars = [x1; x2; x3];
prog = sosprogram(vars);
[prog,V] = sospolyvar(prog,  monomials([x1; x2; x3], 2:4), 'wscoeff');

idx = [10,16,19,4,7,13,9,15,18,1,2,5,11,3,6,12,8,14,17];
K =  [4*coeff_17 3*coeff_18  3*coeff_22  3*coeff_7  2*coeff_19 ...
    2*coeff_23 2*coeff_8   2*coeff_26  2*coeff_11 2*coeff_1 ...
    coeff_20  coeff_24     coeff_9     coeff_27   coeff_12 ...
    coeff_2  coeff_29      coeff_14    coeff_4;
    coeff_18  2*coeff_19   coeff_23    coeff_8   3*coeff_20 ...
    2*coeff_24 2*coeff_9   coeff_27    coeff_12  coeff_2 ...
    4*coeff_21 3*coeff_25  3*coeff_10  2*coeff_28 2*coeff_13 ...
    2*coeff_3  coeff_30    coeff_15   coeff_5;
    coeff_22   coeff_23    2*coeff_26 coeff_11   coeff_24 ...
    2*coeff_27 coeff_12    3*coeff_29 2*coeff_14 coeff_4 ...
    coeff_25   2*coeff_28  coeff_13   3*coeff_30 2*coeff_15 ...
    coeff_5    4*coeff_31  3*coeff_16 2*coeff_6];
K = K(:,idx);
%        K(:,idx)*z - [V1;V2;V3]
%        [V1;V2;V3] - [diff(V,x1); diff(V,x2); diff(V,x3)]
%       V1 - diff(V,x1)

f = F*z + G*K0*z;

prog = sosineq(prog,V);
expr = -[diff(V,x1) diff(V,x2) diff(V,x3)]*f-z'*(Q+K0'*R*K0)*z;%% ;
%V0 =8.518*x1^4 - 4.25e-11*x1^3*x2 - 1.62e-9*x1^3*x3 - 2.971e-10*x1^3 + 16.4*x1^2*x2^2 + 2.018e-9*x1^2*x2*x3 - 1.648e-9*x1^2*x2 + 14.69*x1^2*x3^2 - 9.9e-10*x1^2*x3 + 14.89*x1^2 + 2.051e-9*x1*x2^3 + 3.09e-9*x1*x2^2*x3 + 3.461e-9*x1*x2^2 + 1.344e-9*x1*x2*x3^2 + 5.972*x1*x2*x3 + 1.866e-9*x1*x2 + 1.115e-9*x1*x3^3 + 1.523e-9*x1*x3^2 + 1.766e-9*x1*x3 + 18.78*x2^4 + 2.539e-10*x2^3*x3 + 1.855e-9*x2^3 + 27.06*x2^2*x3^2 + 5.056e-10*x2^2*x3 + 24.83*x2^2 - 1.457e-9*x2*x3^3 + 1.158e-9*x2*x3^2 - 5.191e-10*x2*x3 + 16.76*x3^4 - 4.258e-10*x3^3 + 25.59*x3^2

prog = sosineq(prog,expr);
prog = sossolve(prog);
Vsave = eval(sosgetsol(prog,prog.decvartable));
V0 = sosgetsol(prog,V)
Vold =V0;

% K = sosgetsol(prog,K);
% K = eval(K);
% K = -1/2 * G'*K;
%%
uold = K0*z;


%Ksave(:,:,end+1) = K;
%




vars = [x1; x2; x3];
VEC =  monomials([x1; x2; x3], 2:4);
c = int(int(int(VEC,-1,1),-1,1),-1,1);
for j = 1: 10
    clear prog V
    prog = sosprogram(vars); 
    [prog,V] = sospolyvar(prog, VEC, 'wscoeff');
    K= [4*coeff_17 3*coeff_18  3*coeff_22  3*coeff_7  2*coeff_19 ...
        2*coeff_23 2*coeff_8   2*coeff_26  2*coeff_11 2*coeff_1 ...
        coeff_20  coeff_24     coeff_9     coeff_27   coeff_12 ...
        coeff_2  coeff_29      coeff_14    coeff_4;
        coeff_18  2*coeff_19   coeff_23    coeff_8   3*coeff_20 ...
        2*coeff_24 2*coeff_9   coeff_27    coeff_12  coeff_2 ...
        4*coeff_21 3*coeff_25  3*coeff_10  2*coeff_28 2*coeff_13 ...
        2*coeff_3  coeff_30    coeff_15   coeff_5;
        coeff_22   coeff_23    2*coeff_26 coeff_11   coeff_24 ...
        2*coeff_27 coeff_12    3*coeff_29 2*coeff_14 coeff_4 ...
        coeff_25   2*coeff_28  coeff_13   3*coeff_30 2*coeff_15 ...
        coeff_5    4*coeff_31  3*coeff_16 2*coeff_6];   
    K = K(:, idx);
    prog = sosineq(prog,Vold - V);
    expr = -[diff(V,x1) diff(V,x2) diff(V,x3)]*f-z'*Q*z-uold'*uold;%% ;
    prog = sosineq(prog,expr);

    
    myObj = prog.decvartable * c;
    
    prog = sossetobj(prog, myObj);
    
    prog = sossolve(prog);
    
    V = sosgetsol(prog,V);
    K = sosgetsol(prog,K);
    
    
    K = eval(K);
    K = -1/2 * G'*K;
    Ksave(:,:,end+1) = K;
    
    Vsave(end+1,:) = eval(sosgetsol(prog,prog.decvartable));
    
    Vold =V;
    f = F*z + G*K*z;  
    
    
end


%% Simulate the system dynamics
disp('Simulating the sysm dynamaics ...')
close all
global flag
flag =1;
T = 0.5;
x0 = [2;-1;3];
tsave = [0];
ysave = [x0'];
for i=0:9
    K0 = Ksave(:,:,i+1);
    if i>=6
        flag =0;
    end
    [t, y] = ode45(@RBsys, [i,i+1]*T, ysave(end,:));
    tsave = [tsave; t];
    ysave = [ysave; y];
end
flag = 0;

K0 = Ksave(:,:,1);
[ti, yi] = ode45(@RBsys, [0,i+1]*T,x0);

figure(1)
subplot(311)
plot(ti, yi(:,1), 'r:', tsave,ysave(:,1), 'b', 'LineWidth', 2)
axis([0 5 -0.5 2.5])
legend('With GADP-based controller', 'With initial controller')
xlabel('time (sec)', 'FontSize', 12)
ylabel('x_1', 'FontSize', 12)
annotation(gcf,'textarrow',[0.246428571428571 0.221649029982363],...
    [0.680952380952381 0.572171563660927],'TextEdgeColor','none',...
    'String',{'1st iteration'}, 'FontSize', 12);
annotation(gcf,'textarrow',[0.714285714285714 0.678571428571427],...
    [0.380952380952381 0.26666666666667],'TextEdgeColor','none',...
    'String',{'7th iteration'}, 'FontSize', 12);
annotation(gcf,'textarrow',[0.416071428571429 0.377389770723102],...
    [0.454761904761905 0.312245862884163],'TextEdgeColor','none',...
    'String',{'3rd iteration'}, 'FontSize', 12);
annotation(gcf,'textarrow',[0.272890484739677 0.29745000670855],...
    [0.222222222222222 0.319941160216504],'TextEdgeColor','none',...
    'String',{'2nd iteration'}, 'FontSize', 12) ;
%saveas(gcf,'Ex4_state1.eps', 'psc2')
%saveas(gcf,'Ex4_state1.pdf')


%figure(2)
subplot(312)
plot(ti, yi(:,2), 'r:', tsave,ysave(:,2),'b', 'LineWidth', 2)
axis([0 5 -2 1])
legend('With GADP-based controller', 'With initial controller')
xlabel('time (sec)', 'FontSize', 12)
ylabel('x_2', 'FontSize', 12)
annotation(gcf,'textarrow',[0.24074074074074 0.214506172839506],...
    [0.482269503546099 0.562647754137116],'TextEdgeColor','none',...
    'String',{'1st iteration'}, 'FontSize', 12);
annotation(gcf,'textarrow',[0.319382716049382 0.293148148148147],...
    [0.55919621749409 0.639574468085106],'TextEdgeColor','none',...
    'String',{'2nd iteration'}, 'FontSize', 12);
annotation(gcf,'textarrow',[0.445987654320987 0.370246913580245],...
    [0.546099290780142 0.645579196217494],'TextEdgeColor','none',...
    'String',{'3rd iteration'}, 'FontSize', 12);
annotation(gcf,'textarrow',[0.724444444444444 0.679135802469134],...
    [0.536784741144414 0.644652765699782],'TextEdgeColor','none',...
    'String',{'7th iteration'}, 'FontSize', 12);
%saveas(gcf,'Ex4_state2.eps', 'psc2')
%saveas(gcf,'Ex4_state2.pdf')


subplot(313)
%figure(3)
plot(ti, yi(:,3), 'r:', tsave,ysave(:,3), 'b', 'LineWidth', 2)
axis([0 5 -2 3])
legend('With GADP-based controller', 'With initial controller')
xlabel('time (sec)', 'FontSize', 12)
ylabel('x_3', 'FontSize', 12)
annotation(gcf,'textarrow',[0.219642857142857 0.210934744268077],...
    [0.20952380952381 0.298362039851404],'TextEdgeColor','none',...
    'String',{'1st iteration'}, 'FontSize', 12);
annotation(gcf,'textarrow',[0.296104770453962 0.2875],...
    [0.560317460317463 0.454761904761905],'TextEdgeColor','none',...
    'String',{'2nd iteration'},'FontSize', 12);
annotation(gcf,'textarrow',[0.410714285714285 0.369642857142856],...
    [0.307142857142861 0.421428571428572],'TextEdgeColor','none',...
    'String',{'3rd iteration'}, 'FontSize', 12);
annotation(gcf,'textarrow',[0.719642857142856 0.683928571428569],...
    [0.564285714285715 0.450000000000004],'TextEdgeColor','none',...
    'String',{'7th iteration'},'FontSize', 12);
%saveas(gcf,'Ex4_state3.eps', 'psc2')
%saveas(gcf,'Ex4_state3.pdf')

% Surf V
numm = 40;
xx1 = linspace(-1,1,numm);
xx2 = linspace(-1,1,numm);

[X,Y] = meshgrid(xx1,xx2);
Z = X;

for i = 1:numm
   for j = 1:numm
        x1 = xx1(i);
        x2 = xx2(i);
        x3 = 0;
    z=[ x1^2;
      x1*x2
       x2^2
      x1*x3
      x2*x3
       x3^2
       x1^3
    x1^2*x2
    x1*x2^2
       x2^3
    x1^2*x3
   x1*x2*x3
    x2^2*x3
    x1*x3^2
    x2*x3^2
       x3^3
       x1^4
    x1^3*x2
  x1^2*x2^2
    x1*x2^3
       x2^4
    x1^3*x3
 x1^2*x2*x3
 x1*x2^2*x3
    x2^3*x3
  x1^2*x3^2
 x1*x2*x3^2
  x2^2*x3^2
    x1*x3^3
    x2*x3^3
       x3^4];
 
   Z(i,j) = Vsave(1,:)*z;
   Zn(i,j) = Vsave(end,:)*z;
    end
end
figure(4)
surf(X,Y,Z)
hold on 
surf(X,Y,Zn)
xlabel('x_1', 'FontSize', 12)
ylabel('x_2', 'FontSize', 12)
zlabel('V(x1,x2,0)')
annotation(gcf,'textarrow',[0.687612208258528 0.579892280071813],...
    [0.907433734939759 0.840963855421687],'TextEdgeColor','none',...
    'String',{'V_0(x1,x2,0)'}, 'FontSize', 12);
annotation(gcf,'textarrow',[0.132854578096948 0.143877917414725],...
    [0.183132530120482 0.322168674698801],'TextEdgeColor','none',...
    'String',{'V_7(x1,x2,0)'}, 'FontSize', 12);
%saveas(gcf,'Ex4_cost.eps', 'psc2')
saveas(gcf,'Ex4_cost.pdf')


% %%
% for i = 1:numm
%    for j = 1:numm
%         x1 = 0;
%         x3 = xx2(i);
%         x2 = xx1(i);
%     z=[ x1^2;
%       x1*x2
%        x2^2
%       x1*x3
%       x2*x3
%        x3^2
%        x1^3
%     x1^2*x2
%     x1*x2^2
%        x2^3
%     x1^2*x3
%    x1*x2*x3
%     x2^2*x3
%     x1*x3^2
%     x2*x3^2
%        x3^3
%        x1^4
%     x1^3*x2
%   x1^2*x2^2
%     x1*x2^3
%        x2^4
%     x1^3*x3
%  x1^2*x2*x3
%  x1*x2^2*x3
%     x2^3*x3
%   x1^2*x3^2
%  x1*x2*x3^2
%   x2^2*x3^2
%     x1*x3^3
%     x2*x3^3
%        x3^4];
%  
%    Z(i,j) = Vsave(1,:)*z;
%    Zn(i,j) = Vsave(end,:)*z;
%     end
% end
% figure(5)
% surf(X,Y,Z)
% hold on 
% surf(X,Y,Zn)
% xlabel('x_1', 'FontSize', 12)
% ylabel('x_2', 'FontSize', 12)
% zlabel('V(x1,x2,0)')
% annotation(gcf,'textarrow',[0.687612208258528 0.579892280071813],...
%     [0.907433734939759 0.840963855421687],'TextEdgeColor','none',...
%     'String',{'V_0(x1,x2,0)'}, 'FontSize', 12);
% annotation(gcf,'textarrow',[0.132854578096948 0.143877917414725],...
%     [0.183132530120482 0.322168674698801],'TextEdgeColor','none',...
%     'String',{'V_7(x1,x2,0)'}, 'FontSize', 12);
% 

