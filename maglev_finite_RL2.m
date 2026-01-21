clc 
clear all
close all

% this script disscusses states over time for maglev system based on this
% reaching law   s(k+1) = s(k)−γ1sign[s(k)]min (|s(k)|/(γ1),|s(k)|^β) .
 
%%Where γ1 ∈ R+ and β ∈ (0,1).

t=0:0.1:10; %with step=0.1 and T=10
n=length(t);
z=zeros(3,n); % transformed variable 
x=zeros(3,n); % original variable: position, velocity and current
eq=ones(3,n); % matrix for equilibrium conditions
u=zeros(1,n); % control input
y=zeros(1,n); % measurement quation

phi=[1 0.1 0.005; 0 1 0.1; 0 0 1;];
big_gamma=[0.0002;0.005;0.1];

omega=0.5;
gamma=0.5;
m=0.01187;
Q=0.00018;
beta=0.1;
c=([0.66;1;0.12]);

z(1,1)=0.0155;
z(2,1)=0;
z(3,1)=-18.408;



for k=1:n
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
 plot(t,x(1,:),".--r")
 plot(t,x1d,"--b")
 ylabel('Position(m)')
 xlabel('Time')
 title("Position versus Time")
 axis([0 10 -0.4 0.05]);
 legend('x1(k)','x(1)eq',"Location","best")
 hold off

 % velocity versus time
 subplot(3,2,2);
 x2d=0*ones(1,n);
 hold on
 plot(t,x(2,:),".--r")
 plot(t,x2d,"--b")
 ylabel('Velocity(m/s)')
 xlabel('Time')
 title("Velocity versus Time")
 legend('x2(k)','x(2)eq',"Location","best")
 hold off

 subplot(3,2,3);
 x3d=0.2884*ones(1,n);
 hold on
 plot(t,x(3,:),".--r")
 plot(t,x3d,"--b")
 ylabel('Current(A)')
 xlabel('Time')
 title("Current w.r.t Time")
 legend('x3(k)','x(3d)eq',"Location","best")
 hold off

 subplot(3,2,4);
 x2d=0*ones(1,n);
 hold on
 plot(t,s(1,:),".--r")
 ylabel('sliding variable s(k)')
 xlabel('Time')
 title("sliding variable versus Time")
 ylim([-2 0.5])
 legend('s(k)',"Location","best")
 hold off

 subplot(3,2,5);
 hold on
 plot(t,u(1,:),".--r")
 ylabel('Control Input')
 xlabel('Time')
 title("Control Input versus Time")
 legend('u(t)',"Location","best")
 hold off

