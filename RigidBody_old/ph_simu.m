% Learning interval
global K0
T = 0.25;
x0 = [2;-1;3];
tsave = [0];
ysave = [x0'];
for i=0:9
    K0 = Ksave(:,:,i+1);
    [t, y] = ode45(@RBsys, [i,i+1]*T, ysave(end,:));
    tsave = [tsave; t];
    ysave = [ysave; y];
end

K0 = Ksave(:,:,1);
[ti, yi] = ode45(@RBsys, [0,i+1]*T,x0);

plot(tsave,ysave(:,3), ti, yi(:,3))

