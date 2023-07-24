clear all; close all
%% Constants
lr = 1.5; % in meters
lf = 1.5; % in meters
V = 1; % in m/s
%% Initial condition
% x = [X, Y, PHI]'
xi = [0;0;0];
%% Time limit
tf = 52;
opts = odeset('MaxStep',0.01);
%% A: df is nonzero
df_t = @(t) 45*pi/180;
tspan=0:0.0001:51.78;
[t,x] = ode45(@(t, x) prob4_4ODE(t,x,lf,lr,V,df_t), tspan, xi,opts);
figure;
plot(x(:,1),x(:,2))
xlabel('x (m)'); ylabel('y (m)'); title('Constant df')
grid on
pbaspect([1,1,1])
%% B: df is sinusoid
for V = [1,5]
    for A = [1, 0.5, 0.1]
        df_t = @(t) A*sin(t);
        [t,x] = ode45(@(t, x) prob4_4ODE(t,x,lf,lr,V,df_t), [0, tf], xi,opts);
        figure;
        plot(x(:,1),x(:,2))
        xlabel('x (m)'); ylabel('y (m)'); title(sprintf('Sinusoid df (A=%2.1f, V = %2.1f)',A,V))
        grid on
    end
end
V = 1;
%% C: df is square
amp = 40;
freq = 1/5;
df_t = @(t) -(mod(floor(t*freq),2)*amp-amp/2);
[t,x] = ode45(@(t, x) prob4_4ODE(t,x,lf,lr,V,df_t), [0, tf], xi, opts);
figure; 
plot(x(:,1),x(:,2))
xlabel('x (m)'); ylabel('y (m)'); title('Square df')
grid on

