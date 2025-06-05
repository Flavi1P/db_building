function [pnumsA_matrix, pnumsB_matrix, dataA_matrix, dataB_matrix, ...
          r2s_matrix, ndepths_matrix, depth_matrix] = ...
    RISC_get_full_matchup_regression_data(matchup_pnumsA, matchup_pnumsB, ...
    binned_dataA, binned_dataB, bin_pnumsA, bin_pnumsB, ...
    matchup_r2s, matchup_ndepths, depth_bin_edges)
% RISC_get_full_matchup_regression_data
%
% Extracts vertically-binned profile data for the best (highest r²) matchup
% for each profile in platform A.
%
% INPUTS:
%   matchup_pnumsA     - Profile numbers from platform A [n_matchups x 1]
%   matchup_pnumsB     - Corresponding matched profile numbers from B
%   binned_dataA       - [n_depths x n_profiles] vertically binned data (platform A)
%   binned_dataB       - [n_depths x n_profiles] vertically binned data (platform B)
%   bin_pnumsA         - Profile numbers for columns of binned_dataA
%   bin_pnumsB         - Profile numbers for columns of binned_dataB
%   matchup_r2s        - r² values for each matched pair
%   matchup_ndepths    - Number of valid overlapping depth bins for each pair
%   depth_bin_edges    - Depth bin edges (e.g., 5:5:1000)
%
% OUTPUTS:
%   pnumsA_matrix      - [n_depths x n_profiles] matrix of profile A numbers
%   pnumsB_matrix      - [n_depths x n_profiles] matrix of matched profile B numbers
%   dataA_matrix       - [n_depths x n_profiles] matrix of binned data (A)
%   dataB_matrix       - [n_depths x n_profiles] matrix of binned data (B)
%   r2s_matrix         - [n_depths x n_profiles] matrix of r² values
%   ndepths_matrix     - [n_depths x n_profiles] matrix of n-depths per matchup
%   depth_matrix       - [n_depths x n_profiles] matrix of bin center depths

% ---------------------------------------------
% Step 1: Identify unique profile numbers from platform A
% ---------------------------------------------
unique_pnumsA = unique(matchup_pnumsA(:))';  % row vector of unique profile A numbers
depth_centers = (depth_bin_edges(1:end-1) + depth_bin_edges(2:end)) / 2;
depth_centers = depth_centers(:);            % column vector

n_profiles = length(unique_pnumsA);
n_depths = length(depth_centers);

% ---------------------------------------------
% Step 2: Preallocate output matrices
% ---------------------------------------------
pnumsA_matrix    = repmat(unique_pnumsA, n_depths, 1);         % [n_depths x n_profiles]
depth_matrix     = repmat(depth_centers, 1, n_profiles);       % [n_depths x n_profiles]
pnumsB_matrix    = nan(n_depths, n_profiles);
dataA_matrix     = nan(n_depths, n_profiles);
dataB_matrix     = nan(n_depths, n_profiles);
r2s_matrix       = nan(n_depths, n_profiles);
ndepths_matrix   = nan(n_depths, n_profiles);

% ---------------------------------------------
% Step 3: Loop over each unique profile from platform A
% ---------------------------------------------
for profile_idx = 1:n_profiles
    pnumA = unique_pnumsA(profile_idx);

    % Find all matchups for this profile from A
    match_ixs = find(matchup_pnumsA == pnumA);

    if isempty(match_ixs)
        continue;
    end

    % Select the matchup with the highest r²
    [best_r2, best_ix_rel] = max(matchup_r2s(match_ixs));
    best_ix = match_ixs(best_ix_rel);
    pnumB = matchup_pnumsB(best_ix);

    % Get indices of the binned profiles in the binned data matrices
    idxA = find(bin_pnumsA == pnumA, 1);
    idxB = find(bin_pnumsB == pnumB, 1);

    % Store values into output matrices
    pnumsB_matrix(:, profile_idx)   = pnumB;
    r2s_matrix(:, profile_idx)      = best_r2;
    ndepths_matrix(:, profile_idx)  = matchup_ndepths(best_ix);
    dataA_matrix(:, profile_idx)    = binned_dataA(:, idxA);
    dataB_matrix(:, profile_idx)    = binned_dataB(:, idxB);
end

end
