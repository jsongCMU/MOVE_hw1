function dxdt = prob5ODE(t, x, A, B1, B2, K, phid_t)
% K is 1x4 matrix for computing df using state feedback
% phid_t is a function of time for phid
%% ODE info
% x = [e1, e1', e2, e2']T
% x' = A*x + B1*df + B2*phid
%% Computation
%fprintf('t: %3.3f\n',t);
phid = phid_t(t);
dxdt = (A-B1*K)*x + B2*phid;

