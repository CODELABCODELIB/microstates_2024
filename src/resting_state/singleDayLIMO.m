function [masks, p_vals, models, A, B] = singleDayLIMO(test_values, jids, varargin)
% Enea Ceolini, Leiden University, 26/05/2021
p = inputParser;
addRequired(p, 'test_values');
addRequired(p, 'jids');
addOptional(p, 'FitMethod', 'OLS', @(x) any(validatestring(x, {'OLS', 'IRLS'})));
addOptional(p, 'nBoot', 1000);

parse(p, test_values, jids, varargin{:});
fitMethod = p.Results.FitMethod;
n_boot = p.Results.nBoot;

%% LIMO DAYS
side = size(jids{1}, 1);
n_subs = length(jids);
n_ch = side * side;
n_time = 1;
n_values = size(test_values, 2);

A = zeros(n_subs, n_time, side, side);
B = ones(n_subs, n_values + 1);
B(:, 1:n_values) = test_values;

for i = 1:n_subs
    A(i, 1, :, :) = jids{i};
end

A = log10(reshape(A, n_subs, n_time, n_ch) + 3.1463e-12);

% A = reshape(A, n_subs, n_time, n_ch);
%% NaN guard
A = A(~isnan(B(:, 1)), :, :);
% B = B(~isnan(B(:, 1)), :);
B(isnan(B)) = 0;

%% LIMO

boot_data = permute(A, [3 2 1]); % zeros(n_ch, n_time, n_subs);  % channels, times, individuals
boot_table = limo_create_boot_table(boot_data, n_boot);

M = zeros(n_ch, n_time, n_values);  % channels, times, n_regressor -> F scores
P = zeros(n_ch, n_time, n_values);  % channels, times, n_regressor -> P values
M_full = zeros(n_ch, n_time);  % channels, times -> F scores
P_full = zeros(n_ch, n_time);  % channels, times -> P values

bootM = zeros(n_ch, n_time, n_boot, n_values); % channels, times, nboot, n_regressor
bootP = zeros(n_ch, n_time, n_boot, n_values); % channels, times, nboot, n_regressor
bootM_full = zeros(n_ch, n_time, n_boot); % channels, times, nboot
bootP_full = zeros(n_ch, n_time, n_boot); % channels, times, nboot

W = ones(n_boot, n_subs);
models = cell(n_ch, 1);

for ch = 1:n_ch  % all channels separately
%     multiWaitbar( 'Channels', ch/2500, 'Color', [0.8 0.8 0.1]);
    Y = A(:, :, ch); % (n_subs, n_times);
    X = B; % (n_subs, n_regressors + 1);
    model = limo_glm(Y, X, 0, 0, n_values, fitMethod, 'Time', 0, n_time);
    
    M_full(ch, :) = model.F;
    P_full(ch, :) = model.p;
    
    for m = 1:n_values
        M(ch, :, m) = model.continuous.F(m);
        P(ch, :, m) = model.continuous.p(m);
    end
    models{ch, 1} = model;
    
    model_boot = limo_glm_boot(Y, X, W, 0, 0, n_values, fitMethod, 'Time', boot_table{1, ch});
    
    for j = 1:n_boot
        bootM_full(ch, :, j) = model_boot.F{j};
        bootP_full(ch, :, j) = model_boot.p{j};
    end
    
    for j = 1:n_boot
        for m = 1:n_values
            bootM(ch, :, j, m) = model_boot.continuous.F{j}(m);
            bootP(ch, :, j, m) = model_boot.continuous.p{j}(m);
        end
    end
    
end

% cluster

nM = find_adj_matrix(50, 1);
MCC = 2;
p = 0.05;
% need the single plus the global one
% got rid of a stupid error check: whatever! line 53 of limo_cluster_correction
masks = cell(n_values + 1, 1);
p_vals = cell(n_values + 1, 1);

for m = 1:n_values
    [masks{m}, p_vals{m}] = limo_cluster_correction(M(:, :, m), ...
                                                    P(:, :, m), ...
                                                    bootM(:, :, :, m), ...
                                                    bootP(:, :, :, m), nM, MCC, p);
end

[masks{n_values + 1}, p_vals{n_values + 1}] = limo_cluster_correction(M_full, ...
                                                    P_full, ...
                                                    bootM_full, ...
                                                    bootP_full, nM, MCC, p);

