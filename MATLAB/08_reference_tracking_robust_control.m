%% 08_reference_tracking_robust_control.m
% Bridge crane reference tracking with observer and input constraints

clear; clc; close all;

% Physical parameters
M = 400;    % Cart mass [kg]
m = 600;    % Load mass [kg]
L = 4;      % Cable length [m]
g = 9.8;    % Gravity [m/s^2]

% Linearized system matrices
a32 = (m * g) / M;
a42 = -((M + m) * g) / (M * L);
b2 = 1 / M;
b4 = -1 / (M * L);

A = [0 1 0 0;
     0 0 a32 0;
     0 0 0 1;
     0 0 a42 0];
B = [0; b2; 0; b4];
C = [1 0 L 0];
D = 0; %#ok<NASGU>

% Controllability and observability
fprintf('Controllability rank: %d / %d\n', rank(ctrb(A,B)), size(A,1));
fprintf('Observability rank: %d / %d\n', rank(obsv(A,C)), size(A,1));

% LQR controller
Q = diag([1 0.1 1 0.1]);
R = 10;
[K, ~, ~] = lqr(A, B, Q, R);

% Observer gain
observer_poles = [-2; -2.5; -3; -3.5];
L_observer = place(A', C', observer_poles)';

% Reference tracking gain
N_bar = -inv(C * inv(A - B*K) * B);

% Simulation setup
t_final = 50;
dt = 0.01;
tspan = 0:dt:t_final;
t_change = 15;

x0 = [0; 0; 0; 0];
x_hat0 = [0; 0; 0; 0];
X0 = [x0; x_hat0];

[t_sim, X_sim] = ode45(@(t, X) crane_ode(t, X, A, B, C, K, L_observer, N_bar, t_change), tspan, X0);

x_sim = X_sim(:, 1:4);
x_hat_sim = X_sim(:, 5:8);
y_sim = (C * x_sim')';

y_ref = -5 * ones(size(t_sim));
y_ref(t_sim >= t_change) = 5;

% Control input
u_sim = zeros(length(t_sim), 1);
for i = 1:length(t_sim)
    if t_sim(i) < t_change
        r_t = -5;
    else
        r_t = 5;
    end
    u_sim(i) = -K * x_hat_sim(i,:)' + N_bar * r_t;
    u_sim(i) = max(min(u_sim(i), 500), -500);
end

u_dot_sim = [0; diff(u_sim) ./ diff(t_sim)];
u_dot_sim = max(min(u_dot_sim, 500), -500);

% Plots
figure;
subplot(3,1,1);
plot(t_sim, y_sim, 'b', 'LineWidth', 1.5); hold on;
plot(t_sim, y_ref, 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Load Position y [m]');
title('Load Position vs Reference'); legend('Actual','Reference'); grid on;

subplot(3,1,2);
plot(t_sim, u_sim, 'k', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Control Input u');
title('Control Input'); grid on;

subplot(3,1,3);
plot(t_sim, u_dot_sim, 'm', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Control Input Rate');
title('Control Input Rate'); grid on;

figure;
plot(t_sim, x_sim, 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('States');
title('System States'); legend('x1','x2','phi','phi dot'); grid on;

% Constraint checks
fprintf('Maximum |u(t)|: %.4f\n', max(abs(u_sim)));
fprintf('Maximum |u_dot(t)|: %.4f\n', max(abs(u_dot_sim)));

overshoot = max(abs(y_sim - y_ref)) / 5 * 100;
fprintf('Approximate maximum tracking error percentage: %.2f%%\n', overshoot);

function dXdt = crane_ode(t, X, A, B, C, K, L_obs, N_bar, t_change)
    n = size(A,1);
    x = X(1:n);
    x_hat = X(n+1:end);

    if t < t_change
        r_t = -5;
    else
        r_t = 5;
    end

    y = C * x;
    y_hat = C * x_hat;

    u = -K * x_hat + N_bar * r_t;
    u = max(min(u, 500), -500);

    x_dot = A * x + B * u;
    x_hat_dot = A * x_hat + B * u + L_obs * (y - y_hat);

    dXdt = [x_dot; x_hat_dot];
end
