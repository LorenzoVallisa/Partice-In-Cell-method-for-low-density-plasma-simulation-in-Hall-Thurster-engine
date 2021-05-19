% Matlab script to plot graphs
% uncomment the following and run the commented
% version of contourf and quiver to have lengths
% on the axis instead of the cell number


% n0 = 1e12;
% Te = 1.;
% Qe = 1.602e-19;
% Eps0= 8.854e-12;
% lambda = sqrt(Eps0*Te/(n0*Qe));
% Lx=15*lambda;
% Ly=9*lambda;
% x=0:lambda:Lx;
% y=0:lambda:Ly;

nplots=4;

bigphi=csvread("Phi.csv"); 
phi=zeros(16,10);
figure()
for n=0:nplots-1
    for i=1:16
        phi(i,:)=bigphi((i-1)*10+1+160*n:i*10+160*n);
    end
    subplot(2,2,n+1);
    %contourf(x,y,phi')
    contourf((phi'))
    colorbar
end

bigvel=csvread("NodesVelocity.csv");
vx=phi; vy=phi;
figure()
for n=0:nplots-1
    for i=1:16
        vx(i,:)=bigvel((i-1)*10+1+160*n:i*10+160*n,1);
        vy(i,:)=bigvel((i-1)*10+1+160*n:i*10+160*n,2);
    end
    subplot(2,2,n+1);
    %quiver(x,y,vx',vy')
    quiver(vx',vy')
end
