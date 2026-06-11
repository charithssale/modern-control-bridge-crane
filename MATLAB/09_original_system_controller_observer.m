%% 09_original_system_controller_observer.m
% Robust and stable control design for the original bridge crane system

clear; clc; close all;

% Parameters
U = 4; V = 5; W = 6; Xc = 7; Y = 8; Z = 9;
M = 100 * (5 + Z);
m = 100 * (1 + Y);
L = 5 + Xc;
g = 9.8;

% Linearized system matrices around equilibrium
A = [0 1 0 0;
     0 0 -(m*g)/M 0;
     0 0 0 1;
     0 0 ((M + m)*g)/(M*L) 0];
B = [0; 1/M; 0; -1/(M*L)];
C = [1 0 L 0];
D = 0; %#ok<NASGU>

% LQR controller
Q = diag([1 1 10 10]);
R = 1;
K = lqr(A, B, Q, R);

% Observer design using dual LQR/Kalman-style gain
Qn = diag([0.01 0.01 0.1 0.1]);
Rn = 0.01;
[L_observer, ~, ~] = lqr(A', C', Qn, Rn);
L_observer = L_observer';

% Simulation
dt = 0.01;
T = 10;
time = 0:dt:T;

x0 = [0.1; 0; 0.05; 0];
x_hat0 = [0; 0; 0; 0];
X0 = [x0; x_hat0];

params.A = A;
params.B = B;
params.C = C;
params.K = K;
params.L_observer = L_observer;

[t, X_out] = ode45(@(t, X) crane_dynamics(t, X, params), time, X0);

x = X_out(:,1:4);
cart_position = x(:,1);
cart_velocity = x(:,2);
pendulum_angle = x(:,3);
pendulum_angular_velocity = x(:,4);
output_y = cart_position + L * sin(pendulum_angle);

figure;
subplot(3,2,1); plot(t, cart_position, 'LineWidth', 1.5); title('Cart Position'); grid on;
subplot(3,2,2); plot(t, cart_velocity, 'LineWidth', 1.5); title('Cart Velocity'); grid on;
subplot(3,2,3); plot(t, pendulum_angle, 'LineWidth', 1.5); title('Pendulum Angle'); grid on;
subplot(3,2,4); plot(t, pendulum_angular_velocity, 'LineWidth', 1.5); title('Pendulum Angular Velocity'); grid on;
subplot(3,2,5); plot(t, output_y, 'LineWidth', 1.5); title('Load Position y'); grid on;

function dX = crane_dynamics(~, X, params)
    A = params.A;
    B = params.B;
    C = params.C;
    K = params.K;
    L_observer = params.L_observer;

    x = X(1:4);
    x_hat = X(5:8);

    u = -K * x_hat;
    y = C * x;
    y_hat = C * x_hat;

    dx = A * x + B * u;
    dx_hat = A * x_hat + B * u + L_observer * (y - y_hat);

    dX = [dx; dx_hat];
end
