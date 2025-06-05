function [matchup_ndepths, matchup_r2s, matchup_ms, matchup_m_confints] = ...
    RISC_evaluate_potential_matchups(matchup_pnumsA, matchup_pnumsB, ...
                                     binned_dataA, binned_dataB, ...
                                     bin_pnumsA, bin_pnumsB)
% RISC_evaluate_potential_matchups
%
% Evaluates potential glider profile matchups by computing:
% - Number of overlapping depth bins
% - Type-II regression slope and its confidence interval
% - Coefficient of determination (r²)
%
% INPUTS:
%   matchup_pnumsA   - [n_matchups x 1] Profile numbers from platform A
%   matchup_pnumsB   - [n_matchups x 1] Profile numbers from platform B
%   binned_dataA     - [n_depths x n_profiles] Binned variable data for platform A
%   binned_dataB     - [n_depths x n_profiles] Binned variable data for platform B
%   bin_pnumsA       - [n_profiles x 1] Profile numbers corresponding to columns of binned_dataA
%   bin_pnumsB       - [n_profiles x 1] Profile numbers corresponding to columns of binned_dataB
%
% OUTPUTS:
%   matchup_ndepths     - [n_matchups x 1] Number of valid overlapping depth bins
%   matchup_r2s         - [n_matchups x 1] Coefficient of determination (r²)
%   matchup_ms          - [n_matchups x 1] Type-II regression slopes
%   matchup_m_confints  - [n_matchups x 1] 95% confidence intervals of the slopes

% Initialize
n_matchups = length(matchup_pnumsA);
matchup_ndepths     = nan(n_matchups, 1);
matchup_r2s         = nan(n_matchups, 1);
matchup_ms          = nan(n_matchups, 1);
matchup_m_confints  = nan(n_matchups, 1);

% Loop over all matchup pairs
for i = 1:n_matchups
    % Get indices of the matching profiles in the binned data matrices
    idx_A = find(bin_pnumsA == matchup_pnumsA(i), 1);
    idx_B = find(bin_pnumsB == matchup_pnumsB(i), 1);

    if isempty(idx_A) || isempty(idx_B)
        if isempty(idx_A)
            disp(['Profile ' num2str(matchup_pnumsA(i)) ' not found in Platform A']);
        end
        if isempty(idx_B)
            disp(['Profile ' num2str(matchup_pnumsB(i)) ' not found in Platform B']);
        end
        continue;
    end

    % Extract binned data for the profile pair
    x = binned_dataA(:, idx_A);
    y = binned_dataB(:, idx_B);

    % Identify valid depth bins (non-NaN in both)
    valid_bins = ~isnan(x) & ~isnan(y);
    matchup_ndepths(i) = sum(valid_bins);

    % Only fit regression if enough valid bins
    if matchup_ndepths(i) > 2
        [matchup_ms(i), ~, matchup_m_confints(i), ~, matchup_r2s(i), ~, ~] = ...
            get_linearfit(x(valid_bins), y(valid_bins));
    end
end

end
