%% 04_stabilizability_detectability.m
% Controllability, observability, stabilizability, and detectability analysis

clear; clc;

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

C_matrices = {
    [1 0 L 0], ... % Load position
    [1 0 0 0], ... % Cart position
    [0 0 1 0]  ... % Cable angle
};
outputs = {'Load Position', 'Cart Position', 'Cable Angle'};

eigenvalues_A = eig(A);
unstable_modes = eigenvalues_A(real(eigenvalues_A) >= 0);

for i = 1:length(outputs)
    C = C_matrices{i};
    fprintf('\n=== Analysis for %s ===\n', outputs{i});

    Co = ctrb(A, B);
    Ob = obsv(A, C);

    fprintf('Rank of controllability matrix: %d / %d\n', rank(Co), size(A,1));
    fprintf('Rank of observability matrix: %d / %d\n', rank(Ob), size(A,1));

    is_stabilizable = true;
    for j = 1:length(unstable_modes)
        lambda = unstable_modes(j);
        if rank([lambda * eye(size(A)) - A, B]) ~= size(A,1)
            is_stabilizable = false;
            break;
        end
    end

    is_detectable = true;
    for j = 1:length(unstable_modes)
        lambda = unstable_modes(j);
        if rank([lambda * eye(size(A)) - A; C]) ~= size(A,1)
            is_detectable = false;
            break;
        end
    end

    fprintf('Is system stabilizable? %s\n', string(is_stabilizable));
    fprintf('Is system detectable? %s\n', string(is_detectable));
end
