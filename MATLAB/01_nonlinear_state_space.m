%% 01_nonlinear_state_space.m
% Design and Control of a Bridge Crane System
% Nonlinear state-space model derivation

clear; clc;

% Parameters
U = 4; V = 5; W = 6; Xc = 7; Y = 8; Z = 9;
M = 100 * (5 + Z);      % Cart mass [kg]
m = 100 * (1 + Y);      % Load mass [kg]
L = 5 + Xc;             % Cable length [m]
g = 9.8;                % Gravity [m/s^2]

% Symbolic state variables
syms x1 x2 phi phi_dot u real

% Nonlinear equations
x1_ddot = (u + m * L * (phi_dot^2 * sin(phi) - g * cos(phi))) / (M + m);
phi_ddot = -(x1_ddot * cos(phi) + g * sin(phi)) / L;

% State vector: x = [x1; x2; phi; phi_dot]
x_dot = [x2;
         x1_ddot;
         phi_dot;
         phi_ddot];

% Display equations
disp('Nonlinear State-Space Equations:');
disp('dx1/dt = x2');
disp(['dx2/dt = ', char(x1_ddot)]);
disp('dphi/dt = phi_dot');
disp(['dphi_dot/dt = ', char(phi_ddot)]);
