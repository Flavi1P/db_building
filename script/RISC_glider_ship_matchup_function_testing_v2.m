% RISC_glider_matchup_function_testing_V2.m
% Script to cross-calibrate multiple ocean profiling gliders based on spatial and temporal matchups.

clear;

% Set path to glider data and change directory
RISC_gliders_path = 'C:/Users/flapet/OneDrive - NOC/Documents/IDAPro/lib/db_building/';
glider_data_path = [RISC_gliders_path 'Data/glider/nc_profiled/'];
ship_data_path = [RISC_gliders_path 'data/CTD/DY180/ctd_2db/'];
crosscal_output_path = [RISC_gliders_path 'output/crosscal/'];
if ~exist(crosscal_output_path,'dir')
    mkdir(crosscal_output_path)
end

cd(glider_data_path)

% Load list of NetCDF profile files
glider_files = dir('*_Profiled.nc');
num_gliders = length(glider_files);

% Define NetCDF variable names for gliders
var_time = 'mtime';
var_lat = 'LATITUDE';
var_lon = 'LONGITUDE';
var_profile_num = 'PROFILE_NUMBER';
var_depth = 'PRES_interp';

% Variables of interest for cross-comparison
glider_vars_to_match = {'TEMP','PRAC_SALINITY','BBP700','CHLA','MOLAR_DOXY'};

% Define NetCDF variable names for ship
var_time2 = 'mtime';
var_lat2 = 'latitude';
var_lon2 = 'longitude';
var_profile_num2 = 'ctd_number';
var_depth2 = 'press';

% Variables of interest for cross-comparison
ctd_vars_to_match = {'temp','psal','turbidity','fluor','oxygen'};



%% Load glider NetCDF data
disp('Loading glider data...')
for glider_idx = 1:num_gliders
    % Make sure load_all_nc_variables.m is in the MATLAB path or current directory
    glider_data{glider_idx} = load_all_nc_variables(glider_files(glider_idx).name);
    
    % Convert time to MATLAB datenum format
    glider_data{glider_idx}.mtime = datenum(1970,1,1,0,0,glider_data{glider_idx}.TIME);
    
    % Interpolate pressure across all timestamps
    valid_pressure = ~isnan(glider_data{glider_idx}.PRES);
    glider_data{glider_idx}.PRES_interp = interp1(...
        glider_data{glider_idx}.mtime(valid_pressure), ...
        glider_data{glider_idx}.PRES(valid_pressure), ...
        glider_data{glider_idx}.mtime);
    
    % Store glider name from TRAJECTORY field
    glider_names{glider_idx} = char(glider_data{glider_idx}.TRAJECTORY);
end
%% Load ship NetCDF data
ctd_summary = load_all_nc_variables_transpose('station_summary_dy180_all.nc');
ctd_summary = struct2table(rmfield(ctd_summary,"pad_variable"));
ctd_summary.mtime = datenum(2024,1,1,0,0,ctd_summary.time_start);
ctd_summary.latitude = ctd_summary.lat; ctd_summary.longitude = ctd_summary.lon;
ctd_summary.ctd_number = nan(height(ctd_summary),1);
%%
cd(ship_data_path)
ctd_files = dir('ctd*_2db.nc');
for ctd_idx = 1:length(ctd_files)
    ctd_number = str2double(ctd_files(ctd_idx).name(11:13));
    ctd_summary.ctd_number(ctd_idx) = ctd_number;
    ctd_data = load_all_nc_variables_transpose(ctd_files(ctd_idx).name);
    ctd_data = struct2table(rmfield(ctd_data,"pad_variable"));
    h = height(ctd_data);
    ctd_data.time_total = ctd_summary.time_start(ctd_idx)+ctd_data.time;
    ctd_data.ctd_number(1:h,1) = ctd_number;
    if ctd_idx == 1
        ctd_data_all = ctd_data;
    else
        ctd_data_all = concatonate_tables_robust(ctd_data_all,ctd_data);
    end
end
%%
ctd_data_all.mtime = datenum(2024,1,1,0,0,ctd_data_all.time_total);
ctd_data_all.cp = -log(ctd_data_all.transmittance/100)/0.25; % calculate beam attenuation
ctd_data_all = ctd_data_all(~isnan(ctd_data_all.temp),:);
%% Compute continuous distance to each gliders
disp('Computing ctd-glider distances...')
% for glider_A_idx = 1:num_gliders
    tA = ctd_data_all.(var_time2);
    latA = ctd_data_all.(var_lat2);
    lonA = ctd_data_all.(var_lon2);
    
    % Mask shallow points (<5 m)
    tA(ctd_data_all.(var_depth2) <= 5) = NaN;
    
    figure(20); clf; hold on;
    legend_entries = {};
    
    for glider_B_idx = 1:num_gliders
        tB = glider_data{glider_B_idx}.(var_time);
        latB = glider_data{glider_B_idx}.(var_lat);
        lonB = glider_data{glider_B_idx}.(var_lon);
        tB(glider_data{glider_B_idx}.(var_depth) <= 5) = NaN;

        dist_between = RISC_gliders_continuous_distance(tA, latA, lonA, tB, latB, lonB);
        dist_field_name = ['dist_to_' char(glider_data{glider_B_idx}.TRAJECTORY)];
        ctd_data_all.(dist_field_name) = dist_between;
        
        valid = ~isnan(dist_between + latA + lonA + tA);
        plot(tA(valid), dist_between(valid));
        legend_entries{glider_B_idx} = rem_(char(glider_data{glider_B_idx}.TRAJECTORY));
    end
    
    ylabel(['Distance from ctd (km)']);
    datetick('x');
    xlim_current = get(gca, 'xlim');
    plot(xlim_current, [2 2], 'k'); % threshold line
    ylim([0 50]);
    legend(legend_entries);
    set(gcf,'paperunits','inches','papersize',[12 8],'paperposition',[0 0 12 8])
    print('-dpng',[crosscal_output_path 'cont_dist_ctd.png'])    
% end
%% Extract median values for each profile
disp('Computing median values for each profile...')
for glider_idx = 1:num_gliders
    field_names = fieldnames(glider_data{glider_idx});
    profile_nums = unique(glider_data{glider_idx}.(var_profile_num));
    profile_nums = profile_nums(profile_nums >= 0);
    num_profiles = length(profile_nums);

    profile_summary{glider_idx} = table();

    % Preallocate fields for valid variables
    for var_idx = 1:length(field_names)
        var_name = field_names{var_idx};
        if length(glider_data{glider_idx}.(var_name)) == length(glider_data{glider_idx}.(var_profile_num))
            profile_summary{glider_idx}.(var_name) = nan(num_profiles, 1);
        end
    end

    % Compute per-profile medians
    for p_idx = 1:num_profiles
        profile_mask = glider_data{glider_idx}.(var_profile_num) == profile_nums(p_idx);
        for var_idx = 1:length(field_names)
            var_name = field_names{var_idx};
            if length(glider_data{glider_idx}.(var_name)) == length(glider_data{glider_idx}.(var_profile_num))
                profile_summary{glider_idx}.(var_name)(p_idx) = nanmedian(glider_data{glider_idx}.(var_name)(profile_mask));
            end
        end
    end
end

%% Remove shallow (bad) profiles with median pressure <= 10 m
disp('Filtering shallow profiles...')
for glider_idx = 1:num_gliders
    profile_summary{glider_idx} = profile_summary{glider_idx}(profile_summary{glider_idx}.(var_depth) > 10, :);
end

%% Identify potential matchups by time and distance
disp('Finding potential matchups between ship and gliders...')
max_distance_km = 20;
max_time_diff_days = 2;

for glider_A_idx = 1        
    tA = ctd_summary.(var_time2);
    latA = ctd_summary.(var_lat2);
    lonA = ctd_summary.(var_lon2);
    pnumA = ctd_summary.(var_profile_num2); 

    for glider_B_idx = 1:num_gliders
        % tA = profile_summary{glider_A_idx}.(var_time);
        % latA = profile_summary{glider_A_idx}.(var_lat);
        % lonA = profile_summary{glider_A_idx}.(var_lon);
        % pnumA = profile_summary{glider_A_idx}.(var_profile_num); pnumA(pnumA < 0) = nan;

        tB = profile_summary{glider_B_idx}.(var_time);
        latB = profile_summary{glider_B_idx}.(var_lat);
        lonB = profile_summary{glider_B_idx}.(var_lon);
        pnumB = profile_summary{glider_B_idx}.(var_profile_num); pnumB(pnumB < 0) = nan;

        % Call custom function to get potential matchups
        [pnumsA, pnumsB, dists, dtimes] = RISC_find_potential_matchups(tA, latA, lonA, pnumA, tB, latB, lonB, pnumB, max_distance_km, max_time_diff_days);

        matchup_tables{glider_A_idx, glider_B_idx} = table(pnumsA, pnumsB, dists, dtimes, ...
            'VariableNames', {'matchup_pnumsA', 'matchup_pnumsB', 'matchup_distances', 'matchup_time_diffs'});
        
        num_matchups(glider_A_idx, glider_B_idx) = height(matchup_tables{glider_A_idx, glider_B_idx});
        num_unique_matchups(glider_A_idx, glider_B_idx) = length(unique(pnumsA));
    end
end

%% Create matchup matrices over distance/time bins
disp('Building matchup matrices...')
distance_bins_km = [0.5 1:5 10 15 20];          % km
time_diff_bins_days = [1:3 6 12 24 48] / 24;    % days

clear matchup_matrices
for glider_A_idx = 1%:num_gliders
    for glider_B_idx = 1:num_gliders
        mt = matchup_tables{glider_A_idx, glider_B_idx};
        matchup_matrices{glider_A_idx, glider_B_idx} = RISC_potential_matchup_matrix(...
            mt.matchup_pnumsA, mt.matchup_pnumsB, ...
            mt.matchup_distances, mt.matchup_time_diffs, ...
            distance_bins_km, time_diff_bins_days);
    end
end

%% Plot heatmaps of potential matchups (by distance vs. time)
figure(3); clf;
for glider_A_idx = 1%:num_gliders
    for glider_B_idx = 1:num_gliders
        subplot(1, num_gliders, (glider_A_idx - 1)*num_gliders + glider_B_idx);
        hm = heatmap(num2str(time_diff_bins_days' * 24), num2str(distance_bins_km'), matchup_matrices{glider_A_idx, glider_B_idx});
        xlabel('\Delta time (hours)');
        ylabel('distance (km)');
        title(rem_(['CTD - ' glider_names{glider_B_idx}]));
    end
end
%%
fig_dimensions = [20 5];
set(gcf,'paperunits','inches','papersize',fig_dimensions,'paperposition',[0 0 fig_dimensions])
print('-dpng',[crosscal_output_path 'ship-glider_dist_dt_heatmaps.png'])

%% Bin glider profiles into vertical depth bins
disp('Binning glider profiles vertically...')
depth_bins_m = 5:5:1000;

clear glider_profiles_binned
for glider_idx = 1:num_gliders
    profile_nums = glider_data{glider_idx}.(var_profile_num);
    profile_depths = glider_data{glider_idx}.(var_depth);
    selected_profiles = profile_summary{glider_idx}.(var_profile_num);

    for var_idx = 1:length(glider_vars_to_match)
        var_name = glider_vars_to_match{var_idx};
        profile_values = glider_data{glider_idx}.(var_name);
        [~, glider_profiles_binned{glider_idx}.(var_name)] = ...
            RISC_bin_potential_matchup_profiles(profile_nums, profile_depths, profile_values, selected_profiles, depth_bins_m);
    end
end

%% Bin ctd profiles into vertical depth bins
disp('Binning ship profiles vertically...')
% depth_bins_m = 5:5:1000;

clear ctd_profiles_binned
% for glider_idx = 1:num_gliders
    profile_nums = ctd_data_all.(var_profile_num2);
    profile_depths = ctd_data_all.(var_depth2);
    selected_profiles = ctd_summary.(var_profile_num2);

    for var_idx = 1:length(ctd_vars_to_match)
        var_name = ctd_vars_to_match{var_idx};
        profile_values = ctd_data_all.(var_name);
        [~, ctd_profiles_binned.(var_name)] = ...
            RISC_bin_potential_matchup_profiles(profile_nums, profile_depths, profile_values, selected_profiles, depth_bins_m);
    end
% end

%% Apply stricter filters to matchups (better calibration confidence)
disp('Filtering matchups with stricter distance/time thresholds...')
max_dist_strict_km = 5;
max_time_strict_days = 1;

clear matchup_tables_reduced
for glider_A_idx = 1%:num_gliders
    for glider_B_idx = 1:num_gliders
        mt = matchup_tables{glider_A_idx, glider_B_idx};
        matchup_tables_reduced{glider_A_idx, glider_B_idx} = mt(...
            mt.matchup_distances < max_dist_strict_km & ...
            mt.matchup_time_diffs < max_time_strict_days, :);
    end
end

%% Evaluate each matchup pair for slope, correlation, and depth overlap
disp('Evaluating matchup quality metrics...')
for glider_A_idx = 1%:num_gliders
    for glider_B_idx = 1:num_gliders
        mt = matchup_tables_reduced{glider_A_idx, glider_B_idx};
        pnumsA = mt.matchup_pnumsA;
        pnumsB = mt.matchup_pnumsB;

        bin_pnumsA = ctd_summary.(var_profile_num2);
        bin_pnumsB = profile_summary{glider_B_idx}.(var_profile_num);

        for var_idx = 1:length(glider_vars_to_match)
            var_name = glider_vars_to_match{var_idx};
            var_name2 = ctd_vars_to_match{var_idx};
            bin_dataA = ctd_profiles_binned.(var_name2);
            bin_dataB = glider_profiles_binned{glider_B_idx}.(var_name);

            [n_depths, r2s, slopes, slope_CI] = ...
                RISC_evaluate_potential_matchups(pnumsA, pnumsB, ...
                bin_dataA, bin_dataB, bin_pnumsA, bin_pnumsB);

            % Add results to table
            mt.([var_name '_ndepths']) = n_depths;
            mt.([var_name '_r2s']) = r2s;
            mt.([var_name '_ms']) = slopes;
            mt.([var_name '_m_confints']) = slope_CI;
        end
        matchup_tables_reduced{glider_A_idx, glider_B_idx} = mt;
    end
end

%% Visualize matchups: r² and number of good depth bins
disp('Plotting r² and depth count heatmaps...')

r2_thresholds = [0.99 0.95 0.9 0.85 0.8 0.75 0.7];
ndepth_thresholds = [150 90 45 20 10 5];

clear matchup_r2_matrices matchup_ndepth_matrices

for glider_A_idx = 1%:num_gliders
    for glider_B_idx = 1:num_gliders
        mt = matchup_tables_reduced{glider_A_idx, glider_B_idx};
        pnumsA = mt.matchup_pnumsA;
        pnumsB = mt.matchup_pnumsB;

        matchup_r2_matrices{glider_A_idx, glider_B_idx} = nan(length(glider_vars_to_match), length(r2_thresholds));
        matchup_ndepth_matrices{glider_A_idx, glider_B_idx} = nan(length(glider_vars_to_match), length(ndepth_thresholds));

        for var_idx = 1:length(glider_vars_to_match)
            var_name = glider_vars_to_match{var_idx};
            r2s = mt.([var_name '_r2s']);
            ndepths = mt.([var_name '_ndepths']);

            matchup_r2_matrices{glider_A_idx, glider_B_idx}(var_idx,:) = RISC_matchup_numbers_by_r2(pnumsA, pnumsB, r2s, r2_thresholds);
            matchup_ndepth_matrices{glider_A_idx, glider_B_idx}(var_idx,:) = RISC_matchup_numbers_by_r2(pnumsA, pnumsB, ndepths, ndepth_thresholds);
        end
    end
end

%% Plot r² heatmaps
figure(1); clf;
for glider_A_idx = 1%:num_gliders
    for glider_B_idx = 1:num_gliders
        subplot(1, num_gliders, (glider_A_idx - 1)*num_gliders + glider_B_idx);
        heatmap(num2str(r2_thresholds'), rem_(glider_vars_to_match), matchup_r2_matrices{glider_A_idx, glider_B_idx});
        xlabel('r^2 threshold');
        title(rem_([glider_names{glider_A_idx} ' - ' glider_names{glider_B_idx}]));
    end
end
set(gcf,'paperunits','inches','papersize',fig_dimensions,'paperposition',[0 0 fig_dimensions])
print('-dpng',[crosscal_output_path 'r2_heatmaps_ctd.png'])
%% Plot number of depth bins heatmaps
figure(4); clf;
for glider_A_idx = 1%:num_gliders
    for glider_B_idx = 1:num_gliders
        subplot(1, num_gliders, (glider_A_idx - 1)*num_gliders + glider_B_idx);
        heatmap(num2str(ndepth_thresholds'), rem_(glider_vars_to_match), matchup_ndepth_matrices{glider_A_idx, glider_B_idx});
        xlabel('Min depth bins');
        title(rem_([glider_names{glider_A_idx} ' - ' glider_names{glider_B_idx}]));
    end
end
set(gcf,'paperunits','inches','papersize',fig_dimensions,'paperposition',[0 0 fig_dimensions])
print('-dpng',[crosscal_output_path 'depth_range_heatmaps_ctd.png'])

%% Extract detailed regression data for well-matched profiles
disp('Extracting regression data from selected matchups...')

min_ndepths = 90; % Minimum overlapping depth bins for regression

% Minimum r² thresholds for regression by variable
min_r2.TEMP = 0.99;
min_r2.PRAC_SALINITY = 0.9;
min_r2.BBP700 = 0.95;
min_r2.CHLA = 0.95;
min_r2.MOLAR_DOXY = 0.95;

clear matchup_regression_data_table
for var_idx = 1:length(glider_vars_to_match)
    var_name = glider_vars_to_match{var_idx};
    var_name2 = ctd_vars_to_match{var_idx};
    for glider_A_idx = 1%:num_gliders
        for glider_B_idx = 1:num_gliders
            mt = matchup_tables_reduced{glider_A_idx, glider_B_idx};
            mt = mt(mt.([var_name '_ndepths']) >= min_ndepths, :);

            pnumsA = mt.matchup_pnumsA;
            pnumsB = mt.matchup_pnumsB;
            r2s = mt.([var_name '_r2s']);
            ndepths = mt.([var_name '_ndepths']);

            bin_pnumsA = ctd_summary.(var_profile_num2);
            bin_pnumsB = profile_summary{glider_B_idx}.(var_profile_num);
            dataA = ctd_profiles_binned.(var_name2);
            dataB = glider_profiles_binned{glider_B_idx}.(var_name);

            [pA_mat, pB_mat, dA_mat, dB_mat, r2_mat, nd_mat, z_mat] = ...
                RISC_get_full_matchup_regression_data(pnumsA, pnumsB, dataA, dataB, bin_pnumsA, bin_pnumsB, r2s, ndepths, depth_bins_m);

            % Save to table
            tbl = table(pA_mat(:), pB_mat(:), dA_mat(:), dB_mat(:), r2_mat(:), nd_mat(:), z_mat(:), ...
                'VariableNames', {'ProfileNumA', 'ProfileNumB', [var_name '_A'], [var_name '_B'], 'r2', 'n_depths', 'depth'});

            matchup_regression_data_table{glider_A_idx, glider_B_idx}.(var_name) = tbl;
        end
    end
end
%% Plot regression and extract calibration slopes/offsets
disp('Fitting regression models and computing cross-calibration parameters...')
clear cross_cal_stats
for var_idx = 1:length(glider_vars_to_match)
    var_name = glider_vars_to_match{var_idx};
    var_name2 = ctd_vars_to_match{var_idx};
    figure(200 + var_idx); clf

    for glider_A_idx = 1%:num_gliders
        for glider_B_idx = 1:num_gliders
            reg_table = matchup_regression_data_table{glider_A_idx, glider_B_idx}.(var_name);
            if isempty(reg_table), continue; end

            r2_mask = reg_table.r2 >= min_r2.(var_name);
            x = reg_table.([var_name '_B']);
            y = reg_table.([var_name '_A']);

            subplot(1, num_gliders, (glider_A_idx - 1)*num_gliders + glider_B_idx); hold on
            plot3c(x, y, reg_table.r2, '.');
            plot([min(x) max(x)], [min(x) max(x)], 'k');

            [~, ~, slope, offset, slope_CI, offset_CI, r2] = ...
                plot_linearfit_II(x(r2_mask), y(r2_mask));

            cross_cal_stats.ms.(var_name)(glider_A_idx, glider_B_idx) = slope;
            cross_cal_stats.bs.(var_name)(glider_A_idx, glider_B_idx) = offset;
            cross_cal_stats.mints.(var_name)(glider_A_idx, glider_B_idx) = slope_CI;
            cross_cal_stats.bints.(var_name)(glider_A_idx, glider_B_idx) = offset_CI;
            cross_cal_stats.r2s.(var_name)(glider_A_idx, glider_B_idx) = r2;

            xlabel(rem_([var_name ' ' glider_names{glider_B_idx}]));
            ylabel(rem_([var_name2 ' ctd']));
            clim([0.7 1]);
        end
    end
    cbar_handle = colorbar;
    ylabel(cbar_handle, 'r^2');
    set(gcf,'paperunits','inches','papersize',[6 16],'paperposition',[0 0 6 16])
    print('-dpng',[crosscal_output_path 'ctd_crosscal_regressions' var_name '.png'])    
end

%% Export slope and bias calibration results to CSV
disp('Exporting cross-calibration slope/bias results to CSV...');

for var_idx = 1:length(glider_vars_to_match)
    var_name = glider_vars_to_match{var_idx};
    
    slope_matrix = cross_cal_stats.ms.(var_name);
    bias_matrix  = cross_cal_stats.bs.(var_name);
    slope_ci     = cross_cal_stats.mints.(var_name);
    bias_ci      = cross_cal_stats.bints.(var_name);
    r2_matrix    = cross_cal_stats.r2s.(var_name);

    % Labels
    glider_labels = rem_(glider_names);

    % Initialize arrays
    GliderA = {};
    GliderB = {};
    Slope = [];
    Bias = [];
    Slope_CI = [];
    Bias_CI = [];
    R2 = [];

    % Loop to fill arrays
    for row_idx = 1
        for col_idx = 1:num_gliders
            GliderA{end+1,1}   = 'ctd_DY180';
            GliderB{end+1,1}   = glider_labels{col_idx};
            Slope(end+1,1)     = slope_matrix(row_idx, col_idx);
            Bias(end+1,1)      = bias_matrix(row_idx, col_idx);
            Slope_CI(end+1,1)  = slope_ci(row_idx, col_idx);
            Bias_CI(end+1,1)   = bias_ci(row_idx, col_idx);
            R2(end+1,1)        = r2_matrix(row_idx, col_idx);
        end
    end

    % Assemble final table
    T = table(GliderA, GliderB, Slope, Bias, Slope_CI, Bias_CI, R2);

    % Write to CSV
    csv_filename = fullfile(crosscal_output_path, ['DY180_ctd_glider_calibration_' var_name '.csv']);
    writetable(T, csv_filename);
    disp(['Saved: ' csv_filename]);
end
