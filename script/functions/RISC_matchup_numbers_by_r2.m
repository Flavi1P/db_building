function n_matchups_by_r2 = RISC_matchup_numbers_by_r2(matchup_pnumsA, matchup_pnumsB, matchup_r2s, r2_limits)
% RISC_matchup_numbers_by_r2
%
% For each r² threshold, counts how many unique platform A profiles
% have at least one valid matchup exceeding that threshold.
%
% INPUTS:
%   matchup_pnumsA  - Vector of profile numbers from platform A (one per matchup)
%   matchup_pnumsB  - Vector of profile numbers from platform B (not used here)
%   matchup_r2s     - Vector of r² values for each matchup
%   r2_limits       - Array of r² thresholds to evaluate (e.g., [0.99 0.95 0.9])
%
% OUTPUT:
%   n_matchups_by_r2 - [1 x n_thresholds] Number of unique platform A profiles with matchups ≥ threshold
%

% Number of r² thresholds
n_thresholds = length(r2_limits);
n_matchups_by_r2 = zeros(1, n_thresholds);

% Loop through each r² threshold
for i = 1:n_thresholds
    threshold = r2_limits(i);

    % Logical index of matchups above the threshold
    valid_matchups = matchup_r2s >= threshold;

    if any(valid_matchups)
        % Count unique profile A matchups meeting the r² threshold
        n_matchups_by_r2(i) = numel(unique(matchup_pnumsA(valid_matchups)));
    end
end
end

