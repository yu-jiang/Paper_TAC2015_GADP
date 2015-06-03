clear all
global K
feasibility


x0=[3 -4]; %[3 2]
X=x0;
G=[0;.65];


T=0.005;

syms x1 x2

cs=[ x1*x1;
    x1*x2;
    x2*x2;
    x1*x1*x1;
    x1*x1*x2;
    x1*x2*x2;
    x2*x2*x2;
    x1*x1*x1*x1;
    x1*x1*x1*x2;
    x1*x1*x2*x2;
    x1*x2*x2*x2;
    x2*x2*x2*x2;];

c = double(int(int(cs,0,1),0,1));

clear x1 x2






c=bphi(1,0)+bphi(-1,0)+bphi(0,-1)+bphi(0,1);
psave=[];
ksave=[K];
Xsave=[];
tsave=[];

ji=[0 199 199 199 199 199 199 199 199 199]
jisum=[0 200 400 600 800 1000 1200 1400 1600];% 100+100+100];
for j=1:4
    Phi=[];
    Xi=[];
    Theta=[];
    for i=0:ji(j+1)
        [t,X]=ode45(@jetsysonline,[i*T,(i+1)*T]+jisum(j)*T,[X(end,1:2) zeros(1,35)]);
        Phi=[Phi;X(end,2+1:2+34)];
        Xi=[Xi;X(end,end)];
        Theta=[Theta; bphi(X(end,1),X(end,2))'-bphi(X(1,1),X(1,2))'];
        Xsave=[Xsave; X(:,1:2)];
        tsave=[tsave; t(:)];
    end
    
    cvx_begin sdp
    variable P(5,5) symmetric
    variable L(9,9) symmetric
    p=[P(1,1) P(2,1)+P(1,2) P(2,2) P(1,3)+P(3,1) P(1,4)+P(4,1)+P(2,3)+P(3,2) ...
        P(1,5)+P(5,1)+P(2,4)+P(4,2) P(2,5)+P(5,2) P(3,3) P(3,4)+P(4,3) ...
        P(3,5)+P(5,3)+P(4,4) P(4,5)+P(5,4) P(5,5)]';
    W=[2*p(1) p(2) 3*p(4) 2*p(5) p(6) 4*p(8) 3*p(9) 2*p(10) p(11);
        p(2)   2*p(3) p(5) 2*p(6) 3*p(7) p(9) 2*p(10) 3*p(11) 4*p(12)];
    P<=P0;
    %lk=inv(Phi'*Phi)*Phi'*(-Xi-Theta*p);
    lk=(Phi'*Phi)\(Phi'*(-Xi-Theta*p));
    l=lk(1:25);
    K=lk(26:34)';
    l==[L(1,1);
        L(1,2)+L(2,1);
        L(2,2);
        L(1,3)+L(3,1);
        L(1,4)+L(4,1)+L(2,3)+L(3,2);
        L(1,5)+L(5,1)+L(2,4)+L(4,2);
        L(2,5)+L(5,2);
        L(1,6)+L(6,1)+L(3,3);
        L(1,7)+L(7,1)+L(2,6)+L(6,2)+L(3,4)+L(4,3);
        L(1,8)+L(8,1)+L(2,7)+L(7,2)+L(3,5)+L(5,3)+L(4,4);
        L(1,9)+L(9,1)+L(2,8)+L(8,2)+L(5,4)+L(4,5);
        L(2,9)+L(9,2)+L(5,5);
        L(3,6)+L(6,3);
        L(3,7)+L(7,3)+L(4,6)+L(6,4);
        L(3,8)+L(8,3)+L(4,7)+L(7,4)+L(5,6)+L(6,5);
        L(3,9)+L(9,3)+L(4,8)+L(8,4)+L(5,7)+L(7,5);
        L(4,9)+L(9,4)+L(5,8)+L(8,5);
        L(5,9)+L(9,5);
        L(6,6);
        L(6,7)+L(7,6);
        L(6,8)+L(8,6)+L(7,7);
        L(6,9)+L(9,6)+L(7,8)+L(8,7);
        L(7,9)+L(9,7)+L(8,8);
        L(9,8)+L(8,9);
        L(9,9)];
    L>=0;
    minimize(c'*p)
    cvx_end
    psave=[psave;p(:)'];
    ksave=[ksave;K];
    P0=P;
    if j==1
        P1=P;
    end
end
%%
[t,x]=ode45(@jetsys0,[t(end) 20], Xsave(end,:));
tsave=[tsave;t(:)];
Xsave=[Xsave; x];


%% Draw the figures
close all
K=[K0 0  0     0     0    0      0             0];
[t,x]=ode45(@jetsys0,[0 20], x0);
figure(1)
subplot(211)

plot(tsave,Xsave(:,1),'b-',t,x(:,1),'r-.', 'Linewidth', 2)
axis([0 8 -7 7])
legend('With GADP-based controller', 'With initial controller')
xlabel('time (sec)')
ylabel('x_1')

subplot(212)
plot(tsave,Xsave(:,2),'b-',t,x(:,2),'r-.', 'Linewidth', 2)
axis([0 8 -20 20])
legend('With GADP-based controller', 'With initial controller')
xlabel('time (sec)')
ylabel('x_2')

%for i=1:length(tsave)
%
%end


%subplot(313)
%plot(t,x*K0(:),'b-')%,t,x(:,2),'r:', 'Linewidth', 2)
state_annotations
saveas(gcf,'Ex2_state.eps', 'psc2');
saveas(gcf,'Ex2_state.pdf');
saveas(gcf,'Ex2_state.jpg');



figure(2)
x1=-2:.25:2;
x2=-4:.25:4;
vn=zeros(length(x1),length(x2));
v1=zeros(length(x1),length(x2));
vs=[];
us=[];
un=[];
kn=vn;
k1=v1;
for i=1:length(x1)
    for j=1:length(x2)
        vn(i,j)=phi(x1(i),x2(j))'*P*phi(x1(i),x2(j));
        v1(i,j)=phi(x1(i),x2(j))'*P1*phi(x1(i),x2(j));
        k1(i,j)=ksave(1,:)*sigma(x1(i),x2(j));
        kn(i,j)=ksave(end,:)*sigma(x1(i),x2(j));
    end
end
surf(x1,x2,vn')
hold on
surf(x1,x2,v1')
hold off
xlabel('x_1', 'FontSize', 12)
ylabel('x_2', 'FontSize', 12)
% Create axes
view(gca,[-40.5 14]);
annotation(gcf,'textarrow',[0.216071428571429 0.174535137214669],...
    [0.845238095238095 0.731440045897881],'TextEdgeColor','none','FontSize',12,...
    'String',{'V_0(x1,x2,0)'});
annotation(gcf,'textarrow',[0.132142857142857 0.159949345986154],...
    [0.140476190476191 0.257882960413087],'TextEdgeColor','none','FontSize',12,...
    'String',{'V_4(x1,x2,0)'});
saveas(gcf,'Ex2_cost.eps', 'psc2');
saveas(gcf,'Ex2_cost.pdf');
saveas(gcf,'Ex2_cost.jpg');


figure(3)
surf(x1,x2,kn')
hold on
surf(x1,x2,k1')
hold off
xlabel('x_1')
ylabel('x_2')
saveas(gcf,'Ex2_control_surf.eps', 'psc2');
saveas(gcf,'Ex2_control_surf.pdf') %, 'psc2');