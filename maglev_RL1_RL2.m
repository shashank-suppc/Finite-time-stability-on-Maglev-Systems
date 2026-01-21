clc 
clear all
close all

% this script disscusses states over time for maglev system based on this
% reaching law  s(k +1) = s(k)−sign[s(k)]min(|s(k)|,ω)

t=0:0.1:10; %with step=0.1 and T=10
n=length(t);
z1=zeros(3,n); % transformed variable 
x1=zeros(3,n); % original variable: position, velocity and current
eq=ones(3,n); % matrix for equilibrium conditions
u1=zeros(1,n); % control input
y=zeros(1,n); % measurement quation

z=zeros(3,n); % transformed variable 
x=zeros(3,n); % original variable: position, velocity and current
u=zeros(1,n); % control input

z(1,1)=0.0155;
z(2,1)=0;
z(3,1)=-18.408;


phi=[1 0.1 0.005; 0 1 0.1; 0 0 1;];
big_gamma=[0.0002;0.005;0.1];

omega=0.5;
gamma=0.5;
m=0.01187;
Q=0.00018;
beta=0.1;
c=([0.66;1;0.12]);

z1(1,1)=0.0155;
z1(2,1)=0;
z1(3,1)=-18.408;



for k=1:n
    %position
    x1(1,k)=z1(1,k) + 0.01;
    %x(2) is velocity
    x1(2,k)=z1(2,k);
    %x(3) is current
    x1(3,k)=(x1(1,k).*((9.81-z1(3,k)).*(m./Q)).^0.5);
    s1(k)=transpose(c)*z(:,k);
    u1(1,k)=-inv((transpose(c)*big_gamma))*(transpose(c)*phi*z1(:,k)-s1(k)+sign(s1(k))*min(abs(s1(k)),omega));
    z1(:,k+1)=phi*z1(:,k)+big_gamma*u1(:,k);

    %position
    x(1,k)=z(1,k) + 0.01;
    %x(2) is velocity
    x(2,k)=z(2,k);
    %x(3) is current
    x(3,k)=(x(1,k).*((9.81-z(3,k)).*(m./Q)).^0.5);
    s(k)=transpose(c)*z(:,k);
    u(1,k)=-inv((transpose(c)*big_gamma))*(transpose(c)*phi*z(:,k)-s(k)+gamma.*sign(s(k))*min(abs(s(k))./gamma,abs(s(k)).^beta));
    z(:,k+1)=phi*z(:,k)+big_gamma*u(:,k);
end

 % position versus time
 figure;
 x1d=0.01*ones(1,n);
 subplot(3,2,1);
 hold on
 plot(t,x1(1,:),"*--g")
 plot(t,x(1,:),"o--r")
 plot(t,x1d,"--b")
 ylabel('Position(m)')
 xlabel('Time')
 title("Position versus Time")
 axis([0 10 -0.4 0.05]);
 legend('RL1',"RL2",'x(1)eq',"Location","best")
 hold off

 % velocity versus time
 subplot(3,2,2);
 x2d=0*ones(1,n);
 hold on
 plot(t,x1(2,:),"*--g")
 plot(t,x(2,:),"o--r")
 plot(t,x2d,"--b")
 ylabel('Velocity(m/s)')
 xlabel('Time')
 title("Velocity versus Time")
legend('RL1',"RL2",'x(2)eq',"Location","best")
 hold off

 subplot(3,2,3);
 x3d=0.2884*ones(1,n);
 hold on
  plot(t,x1(3,:),"*--g")
 plot(t,x(3,:),"o--r")
 plot(t,x3d,"--b")
 ylabel('Current(A)')
 xlabel('Time')
 title("Current w.r.t Time")
legend('RL1',"RL2",'x(3)eq',"Location","best")
 hold off

 subplot(3,2,4);
 x2d=0*ones(1,n);
 hold on
 plot(t,s1(1,:),"*--g")
 plot(t,s(1,:),"o--r")
 ylabel('sliding variable s(k)')
 xlabel('Time')
 title("sliding variable versus Time")
 ylim([-2 0.5])
legend('RL1',"RL2","Location","best")
 hold off

 subplot(3,2,5);
 hold on
  plot(t,u1(1,:),"*--g")
 plot(t,u(1,:),"o--r")
 ylabel('Control Input')
 xlabel('Time')
 title("Control Input versus Time")
legend('RL1',"RL2","Location","best")
 hold off

