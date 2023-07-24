% Plot max acceleration as function of lf and h
% a = g * (lf/h * cos(theta)-sin(theta))
%% Constants
g = 9.91; % in m/s^2
m = 100; % in kg
theta = 30; % in deg
step = 0.25; % resolution for quantization
lf = 0.25:step:1.75; % in m
h = 0.5:step:2; % in m
%% Calculation
a = zeros(length(lf), length(h));
% parameterize h
for h_ind=1:length(h)
    cur_h = h(h_ind);
    for lf_ind=1:length(lf)
        cur_lf = lf(lf_ind);
        a(lf_ind,h_ind) = g*(cur_lf/cur_h * cosd(theta) - sind(theta));
    end
end
%% Plot
figure
hold on
for h_ind=1:length(h)
    plot(lf,a(:,h_ind),':o','DisplayName',sprintf('h = %2.2f m',h(h_ind)));
end
xlabel('lf [m]');
ylabel('a [m/s^2]');
lgd = legend('Location','northwest');
title('Max acceleration before tipping vs lf')
