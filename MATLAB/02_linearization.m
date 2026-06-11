%% 02_linearization.m
% Linearization of bridge crane system around phi = 0, phi_dot = 0

clear; clc;

% Parameters
U = 4; V = 5; W = 6; Xc = 7; Y = 8; Z = 9;
M = 100 * (5 + Z);
m = 100 * (1 + Y);
L = 5 + Xc;
g = 9.8;

syms x1 x2 phi phi_dot u real

x1_ddot = (u + m * L * (phi_dot^2 * sin(phi) - g * cos(phi))) / (M + m);
phi_ddot = -(x1_ddot * cos(phi) + g * sin(phi)) / L;

X = [x1; x2; phi; phi_dot];
X_dot = [x2; x1_ddot; phi_dot; phi_ddot];

A_sym = jacobian(X_dot, X);
B_sym = jacobian(X_dot, u);

eq_point = [0; 0; 0; 0; 0];
A = double(subs(A_sym, [X; u], eq_point));
B = double(subs(B_sym, [X; u], eq_point));
C = [1 0 0 0; 0 0 1 0];
D = [0; 0];

disp('Linearized System Matrices:');
disp('A = '); disp(A);
disp('B = '); disp(B);
disp('C = '); disp(C);
disp('D = '); disp(D);
