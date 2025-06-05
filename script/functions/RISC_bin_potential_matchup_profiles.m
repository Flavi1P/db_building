function [depth_matrix, V_matrix] = RISC_bin_potential_matchup_profiles(pnum, depth, V, profile_nums_to_bin, depth_edges)
% RISC_bin_potential_matchup_profiles
%
% Bins vertical profile data (e.g., temperature, salinity) by depth
% for a set of selected profiles. Returns matrices of binned values.
%
% INPUTS:
%   pnum                 - Vector of profile numbers (same length as V and depth)
%   depth                - Vector of depth values (same length as pnum)
%   V                    - Vector of the variable to bin (e.g., TEMP, CHLA)
%   profile_nums_to_bin  - Array of profile numbers to process
%   depth_edges          - Depth bin edges (e.g., 5:5:1000)
%
% OUTPUTS:
%   depth_matrix         - Matrix of binned depth midpoints [n_bins x n_profiles]
%   V_matrix             - Matrix of binned variable medians [n_bins x n_profiles]

% Initialize sizes
num_profiles = length(profile_nums_to_bin);
num_bins = length(depth_edges) - 1;

% Preallocate output matrices
depth_matrix = nan(num_bins, num_profiles);
V_matrix = nan(num_bins, num_profiles);

% Loop over each profile
for profile_idx = 1:num_profiles
    current_profile = profile_nums_to_bin(profile_idx);

    % Identify valid data points for this profile
    valid_points = pnum == current_profile & ~isnan(V) & ~isnan(depth);
    
    % Bin if there are valid values
    if any(valid_points)
        [depth_matrix(:, profile_idx), V_matrix(:, profile_idx)] = ...
            binmedian_edges(depth(valid_points), V(valid_points), depth_edges);
    end
end

end
