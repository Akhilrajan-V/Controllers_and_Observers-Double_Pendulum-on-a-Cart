%% Luenberger Observer for Linear System


%%System Model 
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

B = [0; 1/M; 0; 1/M*l1; 0; 1/M*l2];

Q = [1000 0 0 0 0 0;
     0 1000 0 0 0 0;
     0 0 1000 0 0 0;
     0 0 0 1000 0 0;
     0 0 0 0 100 0;
     0 0 0 0 0 1000;];
 
R = 0.01;


K = lqr(A,B,Q,R);

C1 = [1 0 0 0 0 0];  % C matrix for X 
C2 = [1 0 0 0 0 0; 0 0 0 0 1 0]; % C matrix for X and theta2
C3 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]; % C matrix for X, theta1 and theat2

D = 0; 
     

 %%Initial Condition
 X0 =[0;
     0;
     3.141;
     0;
     0;
     0;
     0;
     0;
     0;
     0;
     0;
     0;];
     
  a21 = zeros(size(A));
  
  c11 = zeros(size(C1));% Observer C matrix for X
  c22 = zeros(size(C2));% Observer C matrix for X and theta2
  c33 = zeros(size(C3));% Observer C matrix for X, theta1 and theat2

  
  pole1 = place(A',C1',[-20 -21 -22 -27 -25 -23;]);

  L1 = pole1';
    
  C_obs = [C1 c11];
   
  A_obs = [(A-B*K) (B*K);
           a21 (A-(L1*C1));];
        
  B_obs = [B;zeros(size(B))];
             
  % Luengerger Observer X   
   sys1 = ss(A_obs, B_obs, C_obs, D);
   
   % Observer of state X
   figure
   step(sys1)
   title('Step Response of State X');
   figure
   initial(sys1,X0)
   title('Response to Initial condition of State X');
   
  % X and theta2 
  pole2 = place(A',C2',[-20 -21 -22 -27 -25 -23;]);
  L2 = pole2';
  C_obs1 = [C2 c22];
   
  A_obs1 = [(A-B*K) (B*K);
           a21 (A-(L2*C2));];
  
   % Luengerger Observer X and theta 2   
   sys1 = ss(A_obs1, B_obs, C_obs1, D);
   
   % Observer of state X and theat 2 
   figure
   step(sys1)
   title('Step Response of States X and theta2');
   
   figure
   initial(sys1,X0)
   title('Response to Initial condition of States X and theta2');
   
 % X, theta1 and theta2 
  pole3 = place(A',C3',[-20 -21 -22 -27 -25 -23;]);
  L3 = pole3';
  C_obs3 = [C3 c33];
   
  A_obs3 = [(A-B*K) (B*K);
           a21 (A-(L3*C3));];
  
   % Luengerger Observer X and theta 2   
   sys1 = ss(A_obs3, B_obs, C_obs3, D);
   
   % Observer of state X and theat 2 
   figure
   step(sys1)
   title('Step Response of States X, theta1 and theta2');
   
   figure
   initial(sys1,X0)
   title('Response to Initial condition of States X, theta1 and theta2');
  
  