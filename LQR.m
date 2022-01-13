%% LQR for Linear System Model


M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.8;



%%System Model 
A = [0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -(g/l1)-(g*m1)/M*l1 0 -(m2*g)/M*l1 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/M*l2 0 -(g/l2)-(g*m2)/M*l2 0];


B = [0; 1/M; 0; 1/M*l1; 0; 1/M*l2];

Q = [1000 0 0 0 0 0;
     0 1000 0 0 0 0;
     0 0 1000 0 0 0;
     0 0 0 1000 0 0;
     0 0 0 0 100 0;
     0 0 0 0 0 1000;];
 
 R = 0.01;
 
 [k, p, e] = lqr(A,B,Q,R);
 
 C = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0;];
     
 D = 0; 

 %%Initial Condition
 X0 =[0;
     0;
     3.141;
     0;
     0;
     0;];
 
 sys = ss((A - B*k), B, C, D);
 
 initial(sys,X0)% Simulation to Initial Condition

 eig(A)% Stability Check
  
 
 
 