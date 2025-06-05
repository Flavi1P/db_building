function matchup_matrix = RISC_potential_matchup_matrix( ...
    matchup_pnumsA, matchup_pnumsB, matchup_distances, matchup_time_diffs, ...
    distance_limits, time_diff_limits)
% RISC_potential_matchup_matrix
%
% Computes a 2D matrix summarizing the number of unique platform A profiles
% with matchups within varying distance and time thresholds.
%
% Each matrix entry (i,j) corresponds to:
%   # of unique profile A matchups with:
%       distance <= distance_limits(i)
%       time_diff <= time_diff_limits(j)
%
% INPUTS:
%   matchup_pnumsA      - Vector of profile numbers from platform A
%   matchup_pnumsB      - Vector of profile numbers from platform B (not used here)
%   matchup_distances   - Distance values for each matchup [km]
%   matchup_time_diffs  - Time difference values for each matchup [days]
%   distance_limits     - Vector of distance thresholds [km]
%   time_diff_limits    - Vector of time difference thresholds [days]
%
% OUTPUT:
%   matchup_matrix      - [m x n] matrix of unique A profiles passing each threshold pair

% Dimensions
m = length(distance_limits);
n = length(time_diff_limits);

% Initialize output matrix
matchup_matrix = zeros(m, n);

% Loop over each combination of distance and time thresholds
for i = 1:m
    for j = 1:n
        within_threshold = (matchup_distances <= distance_limits(i)) & ...
                           (matchup_time_diffs <= time_diff_limits(j));

        if any(within_threshold)
            % Count unique profile A numbers for this bin
            matchup_matrix(i, j) = numel(unique(matchup_pnumsA(within_threshold)));
        end
    end
end
end

