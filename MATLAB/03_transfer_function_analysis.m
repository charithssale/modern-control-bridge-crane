%% 03_transfer_function_analysis.m
% Transfer function and pole-zero analysis

clear; clc; close all;

% Parameters
U = 4; V = 5; W = 6; Xc = 7; Y = 8; Z = 9;
M = 100 * (5 + Z);
m = 100 * (1 + Y);
L = 5 + Xc;
g = 9.8;

% Linearized matrices
b2 = 1 / M;
b3 = -g * (M + m) / (L * M);
b4 = -1 / (L * M);

A = [0 1 0 0;
     0 0 m*g/M 0;
     0 0 0 1;
     0 0 b3 0];
B = [0; b2; 0; b4];

% Output definitions
C1 = [1 0 L 0];   % Load position y = x1 + L*phi
C2 = [1 0 0 0];   % Cart position
C3 = [0 0 1 0];   % Cable angle
D = 0;

sys_y1 = ss(A, B, C1, D);
sys_y2 = ss(A, B, C2, D);
sys_y3 = ss(A, B, C3, D);

G_y1 = tf(sys_y1);
G_y2 = tf(sys_y2);
G_y3 = tf(sys_y3);

disp('Transfer Function for y1 - Load Position:'); disp(G_y1);
disp('Transfer Function for y2 - Cart Position:'); disp(G_y2);
disp('Transfer Function for y3 - Cable Angle:'); disp(G_y3);

figure; pzmap(sys_y1); title('Pole-Zero Map: Load Position'); grid on;
figure; pzmap(sys_y2); title('Pole-Zero Map: Cart Position'); grid on;
figure; pzmap(sys_y3); title('Pole-Zero Map: Cable Angle'); grid on;

poles_y1 = pole(sys_y1);
poles_y2 = pole(sys_y2);
poles_y3 = pole(sys_y3);

disp('Poles for y1:'); disp(poles_y1);
disp('Poles for y2:'); disp(poles_y2);
disp('Poles for y3:'); disp(poles_y3);

disp(['Is y1 asymptotically stable? ', string(all(real(poles_y1) < 0))]);
disp(['Is y2 asymptotically stable? ', string(all(real(poles_y2) < 0))]);
disp(['Is y3 asymptotically stable? ', string(all(real(poles_y3) < 0))]);
