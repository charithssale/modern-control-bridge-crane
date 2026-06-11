%% 05_proportional_controller_analysis.m
% Stability test using a fixed serial proportional controller

clear; clc; close all;

% Parameters
U = 4; V = 5; W = 6; Xc = 7; Y = 8; Z = 9;
M = 100 * (5 + Z);
m = 100 * (1 + Y);
L = 5 + Xc;
g = 9.8;

b2 = 1 / M;
b3 = -g * (M + m) / (L * M);
b4 = -1 / (L * M);

A = [0 1 0 0;
     0 0 m*g/M 0;
     0 0 0 1;
     0 0 b3 0];
B = [0; b2; 0; b4];
C = [1 0 L 0];
D = 0;

sys_open = ss(A, B, C, D);
G_open = tf(sys_open);

disp('Open-loop transfer function:');
disp(G_open);

Kp_values = linspace(-1000, 1000, 500);
stability_results = [];
max_real_poles = zeros(size(Kp_values));

for idx = 1:length(Kp_values)
    Kp = Kp_values(idx);
    G_closed = feedback(G_open, Kp);
    poles_test = pole(G_closed);
    max_real_poles(idx) = max(real(poles_test));

    if all(real(poles_test) < 0)
        stability_results = [stability_results; Kp]; %#ok<AGROW>
    end
end

if isempty(stability_results)
    fprintf('The system cannot be stabilized using a fixed serial proportional controller.\n');
else
    fprintf('Stabilizing Kp range: %.2f to %.2f\n', min(stability_results), max(stability_results));
end

figure;
plot(Kp_values, max_real_poles, 'LineWidth', 1.5);
xlabel('Kp');
ylabel('Maximum Real Part of Closed-Loop Poles');
title('Stability Analysis for Varying Kp');
grid on;
