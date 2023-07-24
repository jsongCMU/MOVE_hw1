function dxdt = prob4_4ODE(t,x,lf,lr,V,df_t)
% Note: df_t is function of time defining df
%% ODE info
% x = [X, Y, PHI]'
% X' = V*cos(PHI) = V*cos(x(3))
% Y' = V*sin(PHI) = V*sin(x(3))
% PHI' = (V/(lf+lr))*tan(df) = (V/(lf+lr))*tan(df)
df = df_t(t);
dxdt = [V*cos(x(3)); V*sin(x(3)); V/(lf+lr)*tan(df)];
