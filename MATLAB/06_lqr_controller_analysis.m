%% 06_lqr_controller_analysis.m
% LQR optimal controller design and linear/nonlinear response analysis

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
B = [0; b2; 0; b4];

Q1 = diag([0.001 0.001 1 0.001]);
Q2 = diag([1 0.001 0.001 0.001]);
Q_set = {Q1, Q2};
p_values = [10^(-2 + U/10), 10, 10^(-5 + U/10), 10^(-10 + U/10)];

optimal_controllers = {};
closed_loop_eigenvalues = {};
labels = {};

for p_idx = 1:length(p_values)
    p = p_values(p_idx);
    for q_idx = 1:length(Q_set)
        Q = Q_set{q_idx};
        R = p;
        [K, ~, E] = lqr(A, B, Q, R);
        optimal_controllers{end+1} = K; %#ok<SAGROW>
        closed_loop_eigenvalues{end+1} = E; %#ok<SAGROW>
        labels{end+1} = sprintf('p=%.2e, Q%d', p, q_idx); %#ok<SAGROW>

        fprintf('\n%s\n', labels{end});
        disp('Optimal Gain K:'); disp(K);
        disp('Closed-Loop Eigenvalues:'); disp(E);
    end
end

% Root locus-like eigenvalue scatter
figure; hold on;
for i = 1:length(closed_loop_eigenvalues)
    eig_vals = closed_loop_eigenvalues{i};
    scatter(real(eig_vals), imag(eig_vals), 60, 'filled');
end
xlabel('Real Part');
ylabel('Imaginary Part');
title('Closed-Loop Eigenvalue Map for LQR Controllers');
legend(labels, 'Location', 'bestoutside');
grid on;

% Linear response
x0_linear = [(6 + 1); 0; 0.035 * (1 + 5); 0];

for i = 1:length(optimal_controllers)
    K = optimal_controllers{i};
    A_cl = A - B * K;
    [t, x] = ode45(@(t, x) A_cl * x, [0 10], x0_linear);

    figure;
    plot(t, x, 'LineWidth', 1.2);
    xlabel('Time [s]'); ylabel('States');
    title(sprintf('Linear Closed-Loop Response - Controller %d', i));
    legend('x1','x2','phi','phi dot'); grid on;

    u = -K * x';
    figure;
    plot(t, u, 'LineWidth', 1.2);
    xlabel('Time [s]'); ylabel('Control Effort u');
    title(sprintf('Control Effort - Controller %d', i)); grid on;
end

% Nonlinear response
x0_nonlinear = [(6 + 1) * m; 0; 0.035 * (1 + 5); 0];

for i = 1:length(optimal_controllers)
    K = optimal_controllers{i};
    nonlinear_model = @(t, x) bridge_crane_nonlinear(t, x, K, M, m, L, g);
    [t, x] = ode45(nonlinear_model, [0 10], x0_nonlinear);

    figure;
    plot(t, x, 'LineWidth', 1.2);
    xlabel('Time [s]'); ylabel('States');
    title(sprintf('Nonlinear Closed-Loop Response - Controller %d', i));
    legend('x1','x2','phi','phi dot'); grid on;

    u = arrayfun(@(idx) -K * x(idx,:)', 1:length(t))';
    figure;
    plot(t, u, 'LineWidth', 1.2);
    xlabel('Time [s]'); ylabel('Control Effort u');
    title(sprintf('Nonlinear Control Effort - Controller %d', i)); grid on;
end

function dxdt = bridge_crane_nonlinear(~, x, K, M, m, L, g)
    u = -K * x;
    x1_ddot = (u + m * L * (x(4)^2 * sin(x(3)) - g * cos(x(3)))) / (M + m);
    phi_ddot = -(x1_ddot * cos(x(3)) + g * sin(x(3))) / L;
    dxdt = [x(2); x1_ddot; x(4); phi_ddot];
end
