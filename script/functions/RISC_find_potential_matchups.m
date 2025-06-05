function [matchup_pnumsA, matchup_pnumsB, matchup_distances, matchup_time_diffs] = ...
    RISC_find_potential_matchups(timeA, latA, lonA, pnumA, ...
                                  timeB, latB, lonB, pnumB, ...
                                  max_distance, max_time_diff)
% RISC_find_potential_matchups
%
% Identifies profile matchups between two glider platforms based on 
% time and spatial proximity (using great-circle distance).
%
% INPUTS:
%   timeA, timeB         - Time vectors (in days since some epoch)
%   latA, lonA           - Latitude/longitude for platform A
%   latB, lonB           - Latitude/longitude for platform B
%   pnumA, pnumB         - Profile number vectors (same length as lat/lon/time)
%   max_distance         - Maximum distance [km] for a valid matchup
%   max_time_diff        - Maximum time difference [days] for a valid matchup
%
% OUTPUTS:
%   matchup_pnumsA       - Profile numbers from platform A (one per matchup)
%   matchup_pnumsB       - Profile numbers from platform B (matched)
%   matchup_distances    - Distance [km] between matched profiles
%   matchup_time_diffs   - Time differences [days] between matched profiles

% Initialize output arrays
matchup_pnumsA = [];
matchup_pnumsB = [];
matchup_distances = [];
matchup_time_diffs = [];

% Unique profile numbers (excluding NaNs)
unique_pnumsA = unique(pnumA(~isnan(pnumA)));
unique_pnumsB = unique(pnumB(~isnan(pnumB)));

nA = length(unique_pnumsA);
nB = length(unique_pnumsB);

% Preallocate summary profile position/time arrays
profile_latA = nan(nA, 1);
profile_lonA = nan(nA, 1);
profile_timeA = nan(nA, 1);

profile_latB = nan(nB, 1);
profile_lonB = nan(nB, 1);
profile_timeB = nan(nB, 1);

% Compute median time/location for each profile (Platform A)
for i = 1:nA
    p = unique_pnumsA(i);
    idx = pnumA == p;
    profile_latA(i) = nanmedian(latA(idx));
    profile_lonA(i) = nanmedian(lonA(idx));
    profile_timeA(i) = nanmedian(timeA(idx));
end

% Compute median time/location for each profile (Platform B)
for j = 1:nB
    p = unique_pnumsB(j);
    idx = pnumB == p;
    profile_latB(j) = nanmedian(latB(idx));
    profile_lonB(j) = nanmedian(lonB(idx));
    profile_timeB(j) = nanmedian(timeB(idx));
end

% Find matchups
for i = 1:nA
    dt = abs(profile_timeB - profile_timeA(i));
    time_match = dt <= max_time_diff;

    if any(time_match)
        dist_km = latlon_dist(profile_latB, profile_lonB, profile_latA(i), profile_lonA(i));
        spatial_match = time_match & dist_km <= max_distance;

        if any(spatial_match)
            n_matches = sum(spatial_match);
            matchup_pnumsA      = [matchup_pnumsA; repmat(unique_pnumsA(i), n_matches, 1)];
            matchup_pnumsB      = [matchup_pnumsB; unique_pnumsB(spatial_match)];
            matchup_distances   = [matchup_distances; dist_km(spatial_match)];
            matchup_time_diffs  = [matchup_time_diffs; dt(spatial_match)];
        end
    end
end
end

%% -------------------------------------------------------
function dist_km = latlon_dist(lat1, lon1, lat2, lon2)
% latlon_dist
% Computes great-circle distance using the haversine formula.
%
% INPUTS:
%   lat1, lon1  - Vectors of lat/lon in degrees
%   lat2, lon2  - Scalars for comparison point
% OUTPUT:
%   dist_km     - Vector of distances [km] from (lat2, lon2) to (lat1, lon1)

R = 6371;  % Earth's radius in kilometers

% Convert degrees to radians
lat1 = deg2rad(lat1); lat2 = deg2rad(lat2);
lon1 = deg2rad(lon1); lon2 = deg2rad(lon2);

dlat = lat2 - lat1;
dlon = lon2 - lon1;

a = sin(dlat / 2).^2 + cos(lat1) .* cos(lat2) .* sin(dlon / 2).^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));

dist_km = R * c;
end
