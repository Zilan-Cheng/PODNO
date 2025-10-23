clear; clc;

% ---- Parameters ----
s = 421;          % grid size
alpha = 2;        % GRF smoothness
tau = 3;          % GRF scale
N = 10000;        % number of samples
f = ones(s, s);   % forcing term

% ---- Preallocate (N, s, s) ----
coeff_all = zeros(N, s, s, 'single');
sol_all   = zeros(N, s, s, 'single');

% ---- Start parallel pool ----
if isempty(gcp('nocreate'))
    parpool('local');   % or parpool('local', 8) if you want to set cores manually
end

fprintf('Generating %d samples in parallel...\n', N);
tic;

% ---- Parallel loop ----
parfor i = 1:N
    % For reproducibility: unique random seed per worker
    rng(i);

    % Generate Gaussian random field
    norm_a = GRF(alpha, tau, s);

    % Threshold coefficient
    thresh_a = zeros(s, s);
    thresh_a(norm_a >= 0) = 12;
    thresh_a(norm_a < 0) = 3;

    % Solve PDE
    thresh_p = solve_gwf(thresh_a, f);

    % Store to (N, s, s)
    coeff_all(i,:,:) = single(thresh_a);
    sol_all(i,:,:)   = single(thresh_p);
end

toc;

% ---- Save dataset ----
save('darcy_threshold_dataset_N_first.mat', 'coeff_all', 'sol_all', '-v7.3');
fprintf('Data saved to darcy_threshold_dataset_N_first.mat\n');
