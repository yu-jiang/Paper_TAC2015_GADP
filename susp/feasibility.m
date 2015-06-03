clear all;
syms x1 x2 x3 x4 real
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m
kn = 0.1*ks;
% State matrices
A = [ 0 1 0 0;
	[-ks -bs ks bs]/mb ; ...
	0 0 0 1;
	[ks bs -ks-kt -bs]/mw];
B = [ 0 0; 0 10000/mb ; 0 0 ; [kt -10000]/mw];
B = B(:,2);
C = [1 0 0 0; 1 0 -1 0; A(2,:)];


Klqr = lqr(A,B,eye(4),1);

f= A*[x1;x2;x3;x4] + [0;-kn*(x1-x3)^3/mb;0;kn*(x1-x3)^3/mw];

q0 =100*x1^2+x2^2+x3^2+x4^2;
rx = 1;%+(x1^2+x2^2+x3^2+x4^2);
vars = [x1;x2;x3;x4];




prog = sosprogram(vars);



% =============================================
% The Lyapunov function V(x):
[prog,V] = sospolyvar(prog,monomials([x1;x2;x3;x4],2:4),'wscoeff');

myObj= int(int(int(int(V,-.5,.5),-10,10),-.5,.5),-10,10);

prog = sosineq(prog,V-0.0001*(x1^2+x2^2+x3^2+x4^2));
expr = -[diff(V,x1) diff(V,x2) diff(V,x3) diff(V,x4)]*rx*f-rx*q0;
prog = sosineq(prog,expr);
prog = sossolve(prog);
V0 = sosgetsol(prog,V)
Vold =V0;
u_prev = [0*x1];
% Iteration
T = .25;
tsave=[0];
%xsave=[-0.2,1,-0.1,-1];
for i=1:10
    i
    %[t,y] = ode23s(@(t,x) mySusp(t,x,u_prev/rx,1), [i-1 i]*T, xsave(end,:));
    %tsave = [tsave;t];
    %xsave=[xsave;y];
	clear prog V

	prog = sosprogram(vars);
	[prog,V] = sospolyvar(prog,monomials([x1;x2;x3;x4],2:4),'wscoeff');
	prog = sosineq(prog,V0-V);
	prog = sosineq(prog, V);
	u_prev = -1/2*B'*[diff(V0,x1) diff(V0,x2) diff(V0,x3) diff(V0,x4)].';
	q =rx*q0+u_prev'*u_prev;
	expr = -[diff(V,x1) diff(V,x2) diff(V,x3) diff(V,x4)]*(rx*f+B*u_prev)-q;
	prog = sosineq(prog,expr);
	prog = sossetobj(prog, myObj);
	prog = sossolve(prog);
	V0 = sosgetsol(prog,V)
end
% % Form a new SOS
%    [t,y] = ode23s(@(t,x) mySusp(t,x,u_prev/rx,0), [i*T 10], xsave(end,:));
%    tsave = [tsave;t];
%    xsave=[xsave;y];
Vnew = V0;

%%
x0 = [0,0,0,10];
%[t1,y1] = ode23s(@(t,x) mySusp(t,x,u_prev/rx,1), [0 10], zeros(4,1));
[t1,y1] = ode23s(@(t,x) mySusp(t,x,100*u_prev,0), [120 123],x0);
%t1= [t1;t11];
%y1 = [y1;y11];
%[t2,y2] = ode23s(@(t,x) mySusp(t,x,u_prev/rx,0), [120 125], x0);
%t=[t;t2];y=[y;y2];
%t1= [t1;t11];
%y1 = [y1;y11];

%[t,y] = ode45(@mySusp0, [0,  15],x0);
[t,y] = ode23s(@mySusp0, [120, 123],x0);
%t=[t;t2];y=[y;y2];
%%
figure(1)
subplot(222); 
plot(t1,y1(:,1),t,y(:,1), 'r:','linewidth',2);
xlabel('time (sec)', 'FontSize',12)
ylabel('x_1','FontSize',12)
legend('With GADP-based controller', 'With no controller')
%axis([0 5 -.5 .5])

subplot(224); 
plot(t1,y1(:,2),t,y(:,2),'r:','linewidth',2);
%axis([0 5 -5 5])
xlabel('time (sec)', 'FontSize',12)
ylabel('x_2','FontSize',12)
legend('With GADP-based controller', 'With no controller')
%subplot(413); plot(tsave,xsave(:,3),t,y(:,3));
%subplot(414); plot(tsave,xsave(:,4),t,y(:,4));
%%
figure(2)
subplot(211); 
plot(t1,y1(:,3),t,y(:,3), 'r:','linewidth',2);
xlabel('time (sec)', 'FontSize',12)
ylabel('x_3','FontSize',12)
legend('With GADP-based controller', 'With no controller')
%axis([0 5 -.5 .5])

subplot(212); 
plot(t1,y1(:,4),t,y(:,4),'r:','linewidth',2);
%axis([0 5 -5 5])
xlabel('time (sec)', 'FontSize',12)
ylabel('x_4','FontSize',12)
legend('With GADP-based controller', 'With no controller')

%%
% % Create textarrow
% annotation(gcf,'textarrow',[0.218120805369128 0.184563758389262],...
%     [0.85401066098081 0.8272921108742],'TextEdgeColor','none',...
%     'String',{'1st iteration'});
% 
% % Create textarrow
% annotation(gcf,'textarrow',[0.23993288590604 0.186107382550335],...
%     [0.168443496801706 0.148443496801705],'TextEdgeColor','none',...
%     'String',{'1st iteration'});
% 
% % Create textarrow
% annotation(gcf,'textarrow',[0.555369127516778 0.518322147651006],...
%     [0.690831556503198 0.741194029850746],'TextEdgeColor','none',...
%     'String',{'Last iteration'});
% 
% % Create textarrow
% annotation(gcf,'textarrow',[0.560268456375839 0.523221476510066],...
%     [0.223070362473348 0.273432835820895],'TextEdgeColor','none',...
%     'String',{'Last iteration'});
%%
% figure(2)
% subplot(211); 
% plot(tsave,xsave(:,3),t,y(:,3), 'r:','linewidth',2);
% xlabel('time (sec)', 'FontSize',12)
% ylabel('x_3','FontSize',12)
% legend('With GADP-based controller', 'With no controller')
% axis([0 5 -.2 .2])
% 
% subplot(212); 
% plot(tsave,xsave(:,4),t,y(:,4),'r:','linewidth',2);
% axis([0 5 -10 15])
% xlabel('time (sec)', 'FontSize',12)
% ylabel('x_4','FontSize',12)
% legend('With GADP-based controller', 'With no controller')
% %subplot(413); plot(tsave,xsave(:,3),t,y(:,3));
% %subplot(414); plot(tsave,xsave(:,4),t,y(:,4));
% 
% 
% % Create textarrow
% annotation(gcf,'textarrow',[0.218120805369128 0.184563758389262],...
%     [0.85401066098081 0.8272921108742],'TextEdgeColor','none',...
%     'String',{'1st iteration'});
% 
% % Create textarrow
% annotation(gcf,'textarrow',[0.23993288590604 0.186107382550335],...
%     [0.168443496801706 0.148443496801705],'TextEdgeColor','none',...
%     'String',{'1st iteration'});
% 
% % Create textarrow
% annotation(gcf,'textarrow',[0.555369127516778 0.518322147651006],...
%     [0.690831556503198 0.741194029850746],'TextEdgeColor','none',...
%     'String',{'Last iteration'});
% 
% % Create textarrow
% annotation(gcf,'textarrow',[0.560268456375839 0.523221476510066],...
%     [0.223070362473348 0.273432835820895],'TextEdgeColor','none',...
%     'String',{'Last iteration'});
% % export_fig Ex3_state -pdf -transparent
%%
figure(3)
xx1=-.4:.04:.4;
xx2=-5:0.5:5;
vn=zeros(length(xx1),length(xx2));
v1=zeros(length(xx1),length(xx2));
vs=[];
us=[];
un=vn;
ulqr = un;
kn=vn;
k1=v1;
x3=0;
x4=0;
for i=1:length(xx1)
    x1 = xx1(i);
    for j=1:length(xx2)
        x2 = xx2(j);
        vn(i,j)=eval(Vnew);
        v1(i,j)=eval(Vold);
       % un(i,j)=eval(u_prev/rx);
%        u1(i,j)=norm(Klqr*[x1;x2;x3;x4]);
%         k1(i,j)=norm([Kold(1,:)*sigma(x1(i),x2(j)),Kold(2,:)*sigma(x1(i),x2(j))]);
%         kn(i,j)=norm([Knew(1,:)*sigma(x1(i),x2(j)),Knew(2,:)*sigma(x1(i),x2(j))]);
    end
end
surf(xx1,xx2,vn')
hold on
surf(xx1,xx2,v1')
hold off
xlabel('x_1', 'FontSize', 12)
ylabel('x_2', 'FontSize', 12)
view(gca,[-30.5 28]);
% Create textarrow
annotation(gcf,'textarrow',[0.210714285714286 0.174535137214669],...
    [0.895238095238095 0.631440045897884],'TextEdgeColor','none','FontSize',12,...
    'String',{'V_0(x1,x2,0,0)'});

% Create textarrow
annotation(gcf,'textarrow',[0.139285714285714 0.186735060271868],...
    [0.183333333333333 0.386454388984516],'TextEdgeColor','none','FontSize',12,...
    'String',{'V_{10}(x1,x2,0,0)'});
% export_fig Ex3_cost -pdf -transparent

%%
xx1=-10:1:10;
xx2=-10:1:10;
vn=zeros(length(xx1),length(xx2));
v1=zeros(length(xx1),length(xx2));
vs=[];
us=[];
un=vn;
ulqr = un;
kn=vn;
k1=v1;
x3=0;
x4=0;
for i=1:length(xx1)
    x1 = xx1(i);
    for j=1:length(xx2)
        x2 = xx2(j);
        un(i,j)=(eval(u_prev/rx));
        u1(i,j)=(Klqr*[x1;x2;x3;x4]);
%         k1(i,j)=norm([Kold(1,:)*sigma(x1(i),x2(j)),Kold(2,:)*sigma(x1(i),x2(j))]);
%         kn(i,j)=norm([Knew(1,:)*sigma(x1(i),x2(j)),Knew(2,:)*sigma(x1(i),x2(j))]);
    end
end
figure(4)
surf(xx1,xx2,un')
hold on
surf(xx1,xx2,u1')
hold off
%saveas(gcf,'Ex2_cost.eps', 'psc2');
%saveas(gcf,'Ex2_cost.pdf');
%saveas(gcf,'Ex2_cost.jpg');
% 
% 
% figure(3)
% surf(x1,x2,kn')
% hold on
% surf(x1,x2,k1')
% hold off
% xlabel('x_1')
% ylabel('x_2')
% %saveas(gcf,'Ex2_control_surf.eps', 'psc2');
% %saveas(gcf,'Ex2_control_surf.pdf') %, 'psc2');
% %%
% figure(4)
% plot(1:8, psave(:,1),'r-o',1:8, psave(:,3), '-*', 1:8, psave(:,5), '-^', 'linewidth',2)
% axis([0.5,8.5,-2, 17])
% legend('p_1 (x_1^2)', 'p_3 (x_2^2)', 'p_5 (x_1^2x_2)')
% xlabel('Iteration', 'FontSize', 12)
% set(gcf,'PaperPositionMode','auto')
% %saveas(gcf,'Ex2_psave.pdf') %, 'psc2');
% %saveas(gcf,'Ex2_psave.eps', 'psc2');
% %export_fig Ex2_psave -pdf -transparent

%%
% figure(1)
% export_fig Ex3_state1 -pdf -transparent
% figure(2)
% export_fig Ex3_state2 -pdf -transparent
% figure(3)
% export_fig Ex3_cost -pdf -transparent
