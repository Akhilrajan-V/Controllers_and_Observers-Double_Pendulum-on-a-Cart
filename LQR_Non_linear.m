%% Non LINEAR LQG CONTROLLER

%%Initial Condition
 X0 =[0;
     0;
     3.141;
     0;
     0;
     0;];
  

time = 0:0.1:350; 

[t1,y1] = ode45(@lqr_nl,time,X0); %using ode45 function for Non linear ODE

figure
plot(t1,y1(:,1))
title('Response of X (Non Linear)')
xlabel('time(Seconds)')
ylabel('Distance(m)')
grid on

figure
plot(t1,y1(:,3))
title('Response of theta1 (Non Linear)')
xlabel('time(Seconds)')
ylabel('radians')
grid on

figure
plot(t1,y1(:,5))
title('Response of theta2 (Non Linear)')
xlabel('time(Seconds)')
ylabel('radians')
grid on

function dydt = lqr_nl(t,y0)
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.8;

A = [0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -(g/l1)-(g*m1)/M*l1 0 -(m2*g)/M*l1 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/M*l2 0 -(g/l2)-(g*m2)/M*l2 0];


B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

Q=[1000 0 0 0 0 0;
   0 1000 0 0 0 0;
   0 0 1000 0 0 0;
   0 0 0 1000 0 0;
   0 0 0 0 100 0;
   0 0 0 0 0 1000];
 
 R = 0.01;
 
 k = lqr(A,B,Q,R);
 
 F=-k*y0; % Input = Gain * Output
 
 dydt=zeros(6,1);
 
 dydt(1) = y0(2); 
 dydt(2)= (F-m1*g*sin(y0(3))*cos(y0(3))-m1*l1*(y0(4)^2)*sin(y0(3))-m2*g*sin(y0(5))*cos(y0(5))-m2*l2*(y0(6)^2)*sin(y0(5)))/(M+m1*((sin(y0(3)))^2)+m2*((sin(y0(5)))^2));
 dydt(3)= y0(4); 
 dydt(4)= (dydt(2)*cos(y0(3))-g*(sin(y0(3))))/l1; 
 dydt(5)= y0(6); 
 dydt(6)= (dydt(2)*cos(y0(5))-g*(sin(y0(5))))/l2; 
end

