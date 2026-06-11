%% 07_kalman_filter_observer.m
% Kalman filter implementation and observer comparison

clear; clc; close all;

% Parameters
U = 4;
M = 100 * (5 + 9);
m = 100 * (1 + 8);
L = 5 + 7;
g = 9.8;

b2 = 1 / M;
b3 = -g * (M + m) / (L * M);
b4 = -1 / (L * M);

A = [0 1 0 0;
     0 0 m*g/M 0;
     0 0 0 1;
     0 0 b3 0];
B = [0; b2; 0; b4]; %#ok<NASGU>
C = [1 0 L 0];
G = eye(4);

v_var = 1e-7;
w_var_values = [10^(-1 + U/10), 1];

kalman_gains = {};

for w_var = w_var_values
    W = w_var * eye(4);
    V = v_var;

    Co_noise = ctrb(A, G * sqrt(W));
    Ob = obsv(A, C);

    fprintf('\nFor w_var = %.2e:\n', w_var);
    fprintf('Controllable with noise input: %s\n', string(rank(Co_noise) == size(A,1)));
    fprintf('Observable: %s\n', string(rank(Ob) == size(A,1)));

    [P, ~, ~] = care(A', C', G * W * G', V);
    L_kalman = P * C' / V;
    A_est = A - L_kalman * C;

    kalman_gains{end+1} = L_kalman; %#ok<SAGROW>

    disp('Kalman Gain:'); disp(L_kalman);
    disp('Estimator Eigenvalues:'); disp(eig(A_est));
end

x0 = [(6 + 1); 0; 0.035 * (1 + 5); 0];
L1 = [1; 0.5; 0.1; 0.01];

observer_set = {kalman_gains{1}, L1};
observer_names = {'Kalman Filter', 'Manual Observer L1'};

for obs_idx = 1:length(observer_set)
    L_obs = observer_set{obs_idx};
    A_est = A - L_obs * C;

    [t, x_est] = ode45(@(t, x) A_est * x, [0 10], x0);

    figure;
    plot(t, x_est, 'LineWidth', 1.2);
    xlabel('Time [s]'); ylabel('Estimated States');
    title(['Estimator Response - ', observer_names{obs_idx}]);
    legend('x1','x2','phi','phi dot'); grid on;

    eig_est = eig(A_est);
    fprintf('\nObserver: %s\n', observer_names{obs_idx});
    fprintf('Stable: %s\n', string(all(real(eig_est) < 0)));
    disp('Eigenvalues:'); disp(eig_est);

    J = trapz(t, sum(x_est.^2, 2));
    fprintf('Quadratic Error: %.6f\n', J);
end
