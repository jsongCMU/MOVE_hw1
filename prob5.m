clear all; close all;
%% Constants
Vx = 30; % m/s
m = 1573; % kg
Iz = 2873; % kg m^2
lf = 1.1; % m
lr = 1.58; % m
Caf = 80000; % N/rad
Car = 80000; % N/rad
%% Matrices
% A
Ac1 = [0; 0; 0; 0];
Ac2 = [1; -2*(Caf+Car)/(m*Vx); 0; -2*(Caf*lf-Car*lr)/(Iz*Vx)];
Ac3 = [0; 2*(Caf+Car)/m; 0; 2*(Caf*lf-Car*lr)/Iz];
Ac4 = [0; 2*(-Caf*lf+Car*lr)/(m*Vx); 1; -2*(Caf*lf^2+Car*lr^2)/(Iz*Vx)];
A = [Ac1, Ac2, Ac3, Ac4];
% B1, B2
B1 = [0; 2*Caf/m; 0; 2*Caf*lf/Iz];
B2 = [0; -2*(Caf*lf-Car*lr)/(m*Vx)-Vx; 0; -2*(Caf*lf^2+Car*lr^2)/(Iz*Vx)];
%% Straight path
% LQR
Q = [15, 0, 0, 0;
    0, 1, 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 25];
R = 1;
K = lqr(A, B1, Q, R);

% Desired position computation
% On y=-5 from t=[0,16.66] -> (1:17)
% On diagonal from t = [16.66, 317.1226] -> (18:317)
% On y=0 from t=[317.1226, inf] -> (318, N)
N = 500;
T = 0.01; % simulation time step
theta = atan(5/90); % slope in rising portion
VT = Vx * T; % distance traveled in one time step
x_des = zeros(N,1);
y_des = zeros(N,1);
phi_des = zeros(N,1);
phid_des = zeros(N,1);
x_des(1) = 0;
y_des(1) = -5;
for i = 2:17
    x_des(i) = x_des(i-1)+VT;
    y_des(i) = -5;
end
for i = 18:317
    x_des(i) = x_des(i-1)+VT*cos(theta);
    y_des(i) = y_des(i-1)+VT*sin(theta);
    phi_des(i) = theta;
end
for i = 318:N
    x_des(i) = x_des(i-1)+VT;
end

% Phi_des'
for i = 2:(N-1)
    phid_des(i) = (phi_des(i+1)-phi_des(i-1))/(2*T);
end
% plot(phid_des)

% Run simulation
tspan = 0:T:(N-1)*T;
phid_t = @(t) phid_des(int64(t/T+1));
xi = [0;0;0;0];
opts = odeset('MaxStep',T);
[t, x] = ode45(@(t,x) prob5ODE(t, x, A, B1, B2, K, phid_t), tspan, xi, opts);

% Compute, display actual position and error
x_act = zeros(N,1);
y_act = zeros(N,1);
for i = 1:N
    x_act(i) = x_des(i)-x(i,1)*sin(phi_des(i)+x(i,3));
    y_act(i) = y_des(i)+x(i,1)*cos(phi_des(i)+x(i,3));
end
figure;
plot(x_des, y_des, 'k:', x_act, y_act, 'b');
title('Deisred and actual vs time (straight)');
figure;
plot(t,x(:,1),'b',t,x(:,3),'r');
title('Error vs time (straight)');

% Compute and display delta and delta'
delta = zeros(N,1);
deltad = zeros(N,1);
for i = 1:N
    delta(i) = -K*x(i,:)';
    if((i>1) && (i < N))
        deltad(i) = (delta(i+1)-delta(i-1))/(2*T);
    end
end
figure;
plot(t,deltad);
title('delta'' vs time (straight)');

% Pass / Fail
e2_max = 0.01;
e1_max_trans = 0.002;
e2_max_trans = 0.0007;
deltad_max = 25;
if(max(abs(deltad))>deltad_max)
    fprintf('FAIL (straight): deltad max too high (%3.3f > %3.3f)\n', max(abs(deltad)),deltad_max)
else
    fprintf('PASS (straight): deltad within limits (%3.3f)\n', max(abs(deltad)))
end
if(max(abs(x(:,3)))>e2_max)
    fprintf('FAIL (straight): e2 max too high (%3.3f > %3.3f)\n', max(abs(x(:,3))),e2_max)
else
    fprintf('PASS (straight): e2 within limits (%3.3f)\n', max(abs(x(:,3))))
end
if(max(abs(x(17,1)-x(117,1)), abs(x(317,1)-x(417,1)))>e1_max_trans)
    fprintf('FAIL (straight): e1 transient too high (%3.3f > %3.3f)\n', max(abs(x(17,1)-x(117,1)), abs(x(317,1)-x(417,1))), e1_max_trans)
else
    fprintf('PASS (straight): e1 transient within limits (%3.5f)\n', max(abs(x(17,1)-x(117,1)), abs(x(317,1)-x(417,1))))
end
if(max(abs(x(17,3)-x(117,3)), abs(x(317,3)-x(417,3)))>e2_max_trans)
    fprintf('FAIL (straight): e2 transient too high (%3.3f > %3.3f)\n\n', max(abs(x(17,3)-x(117,3)), abs(x(317,3)-x(417,3))), e2_max_trans)
else
    fprintf('PASS (straight): e2 transient within limits (%3.5f)\n\n', max(abs(x(17,3)-x(117,3)), abs(x(317,3)-x(417,3))))
end
%% Curved path
% Phid_des = 0 for 1 sec, Vx/1000 for 5 sec, 0 for 1 sec, -Vx/500 for 5 sec
N = 1500;
phid_des = 0:T:(T*(N-1));
for i = 1:100
    phid_des(i) = 0;
end
R = 1000;
for i = 101:600
    phid_des(i) = Vx/R;
end
for i = 601:700
    phid_des(i) = 0;
end
R = -500;
for i = 701:1200
    phid_des(i) = Vx/R;
end
phid_des(1201:end) = 0;
% Phi_des is integral, or running sum, of phid_des
phi_des = 0:T:(T*(N-1));
for i = 2:length(phi_des)
    phi_des(i) = phi_des(i-1) + T*phid_des(i);
end
% Compute x_des and y_des using phi_des and Vx
x_des = 0:T:(T*(N-1));
y_des = 0:T:(T*(N-1));
for i = 2:length(x_des)
    x_des(i) = x_des(i-1) + Vx*cos(phi_des(i));
    y_des(i) = y_des(i-1) + Vx*sin(phi_des(i));
end
% Simulation Setup
Q = [15, 0, 0, 0;
    0, 10, 0, 0;
    0, 0, 15, 0;
    0, 0, 0, 3];
Q = Q;
R = 1;
K = lqr(A, B1, Q, R);
% Simulation
tspan = 0:T:(N-1)*T;
phid_t = @(t) phid_des(int64(t/T+1));
xi = [0;0;0;0];
opts = odeset('MaxStep',T);
[t, x] = ode45(@(t,x) prob5ODE(t, x, A, B1, B2, K, phid_t), tspan, xi, opts);
% Compute actual position
x_act = zeros(N,1);
y_act = zeros(N,1);
for i = 1:N
    x_act(i) = x_des(i)-x(i,1)*sin(phi_des(i)+x(i,3));
    y_act(i) = y_des(i)+x(i,1)*cos(phi_des(i)+x(i,3));
end
% Compute delta and deltad
delta = zeros(N,1);
deltad = zeros(N,1);
for i = 1:N
    delta(i) = -K*x(i,:)';
    if((i>1) && (i < N))
        deltad(i) = (delta(i+1)-delta(i-1))/(2*T);
    end
end
% Display results
figure;
plot(...
    x_des(1:100),y_des(1:100),'k:',...
    x_des(101:600),y_des(101:600),'r:',...
    x_des(601:700),y_des(601:700),'k:',...
    x_des(701:1200),y_des(701:1200),'r:',...
    x_des(1201:end),y_des(1201:end),'k:',...
    x_act, y_act, 'b')
title('Deisred and actual vs time (Curved)');
figure;
plot(t,x(:,1),'b',t,x(:,3),'r');
title(sprintf('Error vs time (Curved), V = %2.2f m/s', Vx));
% Pass / Fail
e1_max = 0.01;
e2_max = 0.01;
deltad_max = 25;
if(max(abs(deltad))>deltad_max)
    fprintf('FAIL (curved): deltad max too high (%3.3f > %3.3f)\n', max(abs(deltad)),deltad_max)
else
    fprintf('PASS (curved): deltad within limits (%3.3f)\n', max(abs(deltad)))
end
if(max(abs(x(:,1)))>e1_max)
    fprintf('FAIL (curved): e1 max too high (%3.3f > %3.3f)\n', max(abs(x(:,1))),e1_max)
else
    fprintf('PASS (curved): e1 within limits (%3.5f)\n', max(abs(x(:,1))))
end

if(max(abs(x(:,3)))>e2_max)
    fprintf('FAIL (curved): e2 max too high (%3.3f > %3.3f)\n', max(abs(x(:,3))),e2_max)
else
    fprintf('PASS (curved): e2 within limits (%3.5f)\n', max(abs(x(:,3))))
end
%% New velocity
Vx = 60; % m/s
% A
Ac1 = [0; 0; 0; 0];
Ac2 = [1; -2*(Caf+Car)/(m*Vx); 0; -2*(Caf*lf-Car*lr)/(Iz*Vx)];
Ac3 = [0; 2*(Caf+Car)/m; 0; 2*(Caf*lf-Car*lr)/Iz];
Ac4 = [0; 2*(-Caf*lf+Car*lr)/(m*Vx); 1; -2*(Caf*lf^2+Car*lr^2)/(Iz*Vx)];
A = [Ac1, Ac2, Ac3, Ac4];
% B1, B2
B1 = [0; 2*Caf/m; 0; 2*Caf*lf/Iz];
B2 = [0; -2*(Caf*lf-Car*lr)/(m*Vx)-Vx; 0; -2*(Caf*lf^2+Car*lr^2)/(Iz*Vx)];
tspan = 0:T:(N-1)*T;
phid_t = @(t) phid_des(int64(t/T+1));
xi = [0;0;0;0];
opts = odeset('MaxStep',T);
[t, x60] = ode45(@(t,x) prob5ODE(t, x, A, B1, B2, K, phid_t), tspan, xi, opts);
figure;
plot(t,x(:,1),'b',t,x(:,3),'r',...
    t,x60(:,1),'b:',t,x60(:,3),'r:');
title('Vx = 30 and Vx = 60 Error vs time');

