% SEITR-SI Epidemic Model 
% Authour by: Awortwe Abraham, Mathematical Sciences,UMaT,Ghana.
% Date 18/06/2024.
% Define pareters
 beta = 0.75; % Transmission probability from IV to SH
 omega = 0.375;  % Transmission probability from IH to SV
 epsilon = 0.25; % Recovery rate
 mu1= 3*10^-5; % natural death rate human
 alpha= 0.0001; % disease induced death rate of human
 bv= 0.05; % birth rate of vector
 bh=4.94*10^-5; % birth rate of human
 mu2= 0.0287; % natural death rate vector
 gamma= 0.04; % Progression rate of exposed human 
 delta= 0.04; % treated rate
 Nh = 130200; % Total population of human
 Nv = 45000; % Total population of vectors
 % Initial conditions
 Sh0 = 90000; % Initial number of susceptible human
 Eh0= 30000;% Inyial number of exposed human
 Ih0 = 10000; % Initial number of infected human
 Th0= (1/50)*Ih0;%initail number of treated human
 Rh0 = 0; % Initial number of recovered human
 Sv0 = 40000; % Initail number of susceptible vector
 Iv0 =5000 ; % Initial number of infected vector
 % Time span tspan = [0 250]; 
% Time period for the simulation 
% Initial state vector y0 = [Sh0 Eh0 Ih0 Th0 Rh0 Sv0 Iv0]; 
% ODE function
% Solve ODE
% Time span
tspan=[0 250];

%Initial state 
y0=[Sh0,Eh0,Ih0,Th0,Rh0,Sv0,Iv0];

%ODE function
f=@(t,y)[ bh - ((beta*y(7)*y(1))/Nv )+ mu1*y(1);((beta*y(7)*y(1)/Nv ))-(gamma + mu1)*y(2);(gamma*y(2))-(alpha + delta +mu1)*y(3);(delta * y(3)) - (mu1+alpha+epsilon)*y(4);(epsilon * y(4)) + mu1*y(5);bv-(omega*y(6)*y(3)/Nh)-mu2*y(6);(omega*y(6)*y(3)/Nh)-mu2*y(7)];
%Check the size of y
disp('Size of yo:');
disp(size(y0));

%Evaluate the ODE function at the initial time and state
t0=tspan(1);%Initial time
fval=f(t0,y0);%Evaluate the ODE function

% Check the size of the vector retuned by ODE function
disp('Size of the vector returned by the ODE function:');
disp(size(fval));

%Solve ODE
options=odeset('RelTol',1e-6,'AbsTol',1e-8);
[t, sol]= ode45(f,tspan,y0,options);
hold on 
%Plot the solution
p=plot(t,sol(:,1),'b',t,sol(:,2),'k',t,sol(:,3),'r',t,sol(:,4),'g',t,sol(:,5),'m',t,sol(:,6),'c',t,sol(:,7),'y');
set(p,'linewidth',2)
legend('S_H','E_H','I_H','T_H','R_H','S_V','I_V')
xlabel('Time(days)')
ylabel('Population')