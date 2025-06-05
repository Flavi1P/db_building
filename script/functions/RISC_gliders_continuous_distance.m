function distA = RISC_gliders_continuous_distance(timeA, latA, lonA, timeB, latB, lonB)
% RISC_gliders_continuous_distance
%
% Computes the time-continuous horizontal (great-circle) distance between
% platform A and platform B using interpolation of B's track onto A's time points.
%
% INPUTS:
%   timeA, latA, lonA - Vectors of timestamps (days), latitudes, and longitudes for platform A
%   timeB, latB, lonB - Vectors of timestamps (days), latitudes, and longitudes for platform B
%
% OUTPUT:
%   distA             - Vector of distances [km] at each A timepoint; NaN where B data is not available
%

% Initialize output with NaNs
distA = nan(size(timeA));

% Clean B’s data
validB = ~isnan(timeB) & ~isnan(latB) & ~isnan(lonB);

% If B has no valid data, exit
if all(~validB)
    return;
end

% Filter valid B data
timeB = timeB(validB);
latB = latB(validB);
lonB = lonB(validB);

% Identify overlapping time range
validA = timeA >= min(timeB) & timeA <= max(timeB) & ~isnan(timeA);

% Interpolate B’s lat/lon to A’s time
latB_interp = interp1(timeB, latB, timeA, 'linear', NaN);
lonB_interp = interp1(timeB, lonB, timeA, 'linear', NaN);

% Compute distance only for overlapping timepoints
distA(validA) = latlon_dist(latA(validA), lonA(validA), latB_interp(validA), lonB_interp(validA));
end

%% -------------------------------------------------------
function dist_km = latlon_dist(lat1, lon1, lat2, lon2)
% latlon_dist
% Computes great-circle distance between two lat/lon coordinate sets using the haversine formula.
%
% INPUTS:
%   lat1, lon1 - First point(s) in decimal degrees (can be vectors)
%   lat2, lon2 - Second point(s), same size as lat1/lon1 or scalars
%
% OUTPUT:
%   dist_km    - Distance in kilometers

R = 6371;  % Earth radius in km

% Convert to radians
lat1 = deg2rad(lat1); lat2 = deg2rad(lat2);
lon1 = deg2rad(lon1); lon2 = deg2rad(lon2);

dlat = lat2 - lat1;
dlon = lon2 - lon1;

a = sin(dlat/2).^2 + cos(lat1) .* cos(lat2) .* sin(dlon/2).^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));

dist_km = R * c;
end
