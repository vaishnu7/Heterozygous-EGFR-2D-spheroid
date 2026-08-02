% Run this in MATLAB after your simulation completes

clear all; close all; clc;

%% =========================================================================
% STEP 1: Load all data files
% =========================================================================

data_dir = 'path-to-your\cellData\';

% Check if directory exists
if ~exist(data_dir, 'dir')
    error('Data directory not found: %s\n', data_dir);
end

% Find all data files
file_pattern = [data_dir 'cell_data_time_*.dat'];
file_list = dir(file_pattern);

if length(file_list) == 0
    error('No data files found in %s\n', data_dir);
end

fprintf('Found %d data files\n\n', length(file_list));

% Load all files
all_data = {};
timesteps = [];

for i = 1:length(file_list)
    filename = fullfile(data_dir, file_list(i).name);
    
    % Read the file, skipping the header (first line)
    opts = delimitedTextImportOptions("NumVariables", 17);
    opts.Delimiter = "\t";
    opts.DataLines = [2, Inf];  % Skip header line
    opts.VariableNames = {'Time', 'TimeStep', 'CellIndex', 'DistanceFromCenter', 'PositionX', 'PositionY', ...
                          'CellAge', 'BirthTime', 'Oxygen', 'CellState', 'TimeHypoxic', 'HIF1Alpha', ...
                          'TGFAlpha', 'EGFRActivation', 'DistanceToBoundary', 'MutationState', 'ProliferativeType'};
    
    % Set variable types
    for j = 1:15
        opts.VariableTypes{j} = 'double';
    end
    opts.VariableTypes{16} = 'string';  % MutationState
    opts.VariableTypes{17} = 'string';  % ProliferativeType
    
    data = readtable(filename, opts);
    all_data{i} = data;
    
    if height(data) > 0
        timesteps(i) = data.Time(1);
    end
end

fprintf('Loaded %d timesteps\n', length(all_data));
fprintf('Time range: %.1f to %.1f hours\n\n', min(timesteps), max(timesteps));

%% =========================================================================
% STEP 2: Extract temporal dynamics and SORT BY TIME
% =========================================================================

fprintf('Processing temporal data...\n');

time_data = [];

for i = 1:length(all_data)
    data = all_data{i};
    
    if height(data) > 0
        time_data.time(i) = data.Time(1);
        time_data.total_cells(i) = height(data);
        time_data.viable(i) = sum(data.CellState == 1);
        time_data.hypoxic(i) = sum(data.CellState == 2);
        time_data.necrotic(i) = sum(data.CellState == 3);
        time_data.mean_oxygen(i) = mean(data.Oxygen);
        time_data.mean_hif1a(i) = mean(data.HIF1Alpha);
        time_data.mean_tgfa(i) = mean(data.TGFAlpha);
        time_data.mean_egfr(i) = mean(data.EGFRActivation);
        time_data.mean_age(i) = mean(data.CellAge);
        time_data.spheroid_radius(i) = max(data.DistanceFromCenter);
    end
end

% SORT ALL DATA BY TIME (fixes the scrambled lines issue)
[time_sorted, sort_idx] = sort(time_data.time);

% Re-sort all_data to match time_data
all_data_sorted = all_data(sort_idx);
all_data = all_data_sorted;

time_data.time = time_sorted;
time_data.total_cells = time_data.total_cells(sort_idx);
time_data.viable = time_data.viable(sort_idx);
time_data.hypoxic = time_data.hypoxic(sort_idx);
time_data.necrotic = time_data.necrotic(sort_idx);
time_data.mean_oxygen = time_data.mean_oxygen(sort_idx);
time_data.mean_hif1a = time_data.mean_hif1a(sort_idx);
time_data.mean_tgfa = time_data.mean_tgfa(sort_idx);
time_data.mean_egfr = time_data.mean_egfr(sort_idx);
time_data.mean_age = time_data.mean_age(sort_idx);
time_data.spheroid_radius = time_data.spheroid_radius(sort_idx);

% Calculate growth rate (from sorted data)
growth_rate = diff(time_data.total_cells) ./ diff(time_data.time);

% ===== DIAGNOSTIC: Check TGF-a data =====
fprintf('\n=== TGF-a Diagnostic ===\n');
fprintf('TGF-a raw values (first 10): %s\n', mat2str(time_data.mean_tgfa(1:min(10, length(time_data.mean_tgfa))), 4));
fprintf('TGF-a min: %.6f\n', min(time_data.mean_tgfa));
fprintf('TGF-a max: %.6f\n', max(time_data.mean_tgfa));
fprintf('TGF-a nonzero count: %d / %d\n', sum(time_data.mean_tgfa > 0), length(time_data.mean_tgfa));
fprintf('TGF-a column index in data: 13\n');

% Check a sample file for TGF-a values
if length(all_data) > 0
    sample_data = all_data{end};
    fprintf('\nSample from last timestep:\n');
    fprintf('  TGFAlpha values (first 5): %s\n', mat2str(sample_data.TGFAlpha(1:min(5, height(sample_data))), 4));
    fprintf('  TGFAlpha min: %.6f, max: %.6f\n', min(sample_data.TGFAlpha), max(sample_data.TGFAlpha));
end
fprintf('=====================================\n\n');

fprintf('Population statistics:\n');
fprintf('  Initial cells: %d\n', time_data.total_cells(1));
fprintf('  Final cells: %d\n', time_data.total_cells(end));
fprintf('  Growth: %.1f%%\n\n', 100*(time_data.total_cells(end) - time_data.total_cells(1))/time_data.total_cells(1));

%% =========================================================================
% FIGURE 1: Population Dynamics
% =========================================================================

fprintf('Creating Figure 1: Population Dynamics...\n');

fig1 = figure('Name', 'Population Dynamics', 'Position', [100, 100, 1400, 600]);

% Plot 1: Total and by state
subplot(2, 3, 1);
plot(time_data.time, time_data.total_cells, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Total');
hold on;
plot(time_data.time, time_data.viable, 'g-', 'LineWidth', 2, 'DisplayName', 'Viable');
plot(time_data.time, time_data.hypoxic, 'y-', 'LineWidth', 2, 'DisplayName', 'Hypoxic');
plot(time_data.time, time_data.necrotic, 'r-', 'LineWidth', 2, 'DisplayName', 'Necrotic');
xlabel('Time (hours)', 'FontSize', 11);
ylabel('Cell Count', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 10);
title('Cell Population Dynamics', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
hold off;

% Plot 2: Cell state fractions
subplot(2, 3, 2);
viable_frac = time_data.viable ./ time_data.total_cells * 100;
hypoxic_frac = time_data.hypoxic ./ time_data.total_cells * 100;
necrotic_frac = time_data.necrotic ./ time_data.total_cells * 100;

plot(time_data.time, viable_frac, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Viable');
hold on;
plot(time_data.time, hypoxic_frac, 'y-', 'LineWidth', 2.5, 'DisplayName', 'Hypoxic');
plot(time_data.time, necrotic_frac, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Necrotic');
xlabel('Time (hours)', 'FontSize', 11);
ylabel('Percentage (%)', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 10);
ylim([0, 100]);
title('Cell State Fractions', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
hold off;

% Plot 3: Mean oxygen
subplot(2, 3, 3);
plot(time_data.time, time_data.mean_oxygen, 'b-', 'LineWidth', 2.5);
xlabel('Time (hours)', 'FontSize', 11);
ylabel('Mean Oxygen Concentration', 'FontSize', 11);
title('Mean Oxygen Level', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 4: HIF-1a and EGFR (normalized to 0-1)
subplot(2, 3, 4);

% Normalize HIF-1a and EGFR to 0-1 range
hif1a_norm = time_data.mean_hif1a / max(time_data.mean_hif1a + eps);
egfr_norm = time_data.mean_egfr / max(time_data.mean_egfr + eps);

plot(time_data.time, hif1a_norm, 'r-', 'LineWidth', 2.5, 'DisplayName', 'HIF-1a');
hold on;
plot(time_data.time, egfr_norm, 'g-', 'LineWidth', 2.5, 'DisplayName', 'EGFR');
xlabel('Time (hours)', 'FontSize', 11);
ylabel('Normalized Concentration (0-1)', 'FontSize', 11);
ylim([0, 1.1]);
legend('Location', 'best', 'FontSize', 10);
title('HIF-1a & EGFR (Normalized)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
hold off;

% Plot 5: TGF-a (separate due to small values)
subplot(2, 3, 5);
tgfa_norm = time_data.mean_tgfa / max(time_data.mean_tgfa + eps);
plot(time_data.time, tgfa_norm, 'color', [1 0.5 0], 'LineWidth', 2.5, 'DisplayName', 'TGF-a');
xlabel('Time (hours)', 'FontSize', 11);
ylabel('Normalized Concentration (0-1)', 'FontSize', 11);
ylim([0, 1.1]);
legend('Location', 'best', 'FontSize', 10);
title('TGF-a (Normalized)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 6: Spheroid growth
subplot(2, 3, 5);
plot(time_data.time, time_data.spheroid_radius, 'b-', 'LineWidth', 2.5);
xlabel('Time (hours)', 'FontSize', 11);
ylabel('Radius (cell units)', 'FontSize', 11);
title('Spheroid Growth', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 6: Population growth (smooth)
subplot(2, 3, 6);
window_size = max(3, round(length(time_data.time) / 20));  % Adaptive smoothing
growth_rate_smooth = movmean(growth_rate, window_size);
plot(time_data.time(1:end-1), growth_rate_smooth, 'b-', 'LineWidth', 2.5);
hold on;
fill([time_data.time(1:end-1), fliplr(time_data.time(1:end-1))], ...
     [growth_rate_smooth + movstd(growth_rate, window_size), ...
      fliplr(growth_rate_smooth - movstd(growth_rate, window_size))], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (hours)', 'FontSize', 11);
ylabel('Growth Rate (cells/hour)', 'FontSize', 11);
title('Population Growth Rate', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
hold off;

sgtitle('PC9 Spheroid - Temporal Dynamics', 'FontSize', 14, 'FontWeight', 'bold');

%% =========================================================================
% FIGURE 1B: TGF-a Detailed Analysis
% =========================================================================

% fprintf('Creating Figure 1B: TGF-a Analysis...\n');
% fig1b = figure('Name', 'TGF-a Analysis', 'Position', [100, 100, 800, 600]);
% 
% plot(time_data.time, time_data.mean_tgfa, 'color', [1 0.5 0], 'LineWidth', 2.5);
% xlabel('Time (hours)', 'FontSize', 12);
% ylabel('Mean TGF-a Concentration', 'FontSize', 12);
% title('TGF-a Over Time', 'FontSize', 14, 'FontWeight', 'bold');
% grid on;

fprintf('Creating Figure 1B: TGF-a Analysis...\n');
fig1b = figure('Name', 'TGF-a Analysis', 'Position', [100, 100, 1400, 500]);

% Plot 1: TGF-a over time (raw values)
subplot(1, 3, 1);
plot(time_data.time, time_data.mean_tgfa, 'color', [1 0.5 0], 'LineWidth', 2.5);
xlabel('Time (hours)', 'FontSize', 11);
ylabel('Mean TGF-a Concentration', 'FontSize', 11);
title('TGF-a Over Time (Raw Values)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 2: TGF-a normalized (0-1)
subplot(1, 3, 2);
tgfa_norm = time_data.mean_tgfa / max(time_data.mean_tgfa + eps);
plot(time_data.time, tgfa_norm, 'color', [1 0.5 0], 'LineWidth', 2.5);
xlabel('Time (hours)', 'FontSize', 11);
ylabel('Normalized TGF-a (0-1)', 'FontSize', 11);
ylim([0, 1.1]);
title('TGF-a Over Time (Normalized)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 3: TGF-a vs HIF-1a comparison
subplot(1, 3, 3);
hif1a_norm = time_data.mean_hif1a / max(time_data.mean_hif1a + eps);
tgfa_norm = time_data.mean_tgfa / max(time_data.mean_tgfa + eps);
plot(time_data.time, hif1a_norm, 'r-', 'LineWidth', 2.5, 'DisplayName', 'HIF-1a');
hold on;
plot(time_data.time, tgfa_norm, 'color', [1 0.5 0], 'LineWidth', 2.5, 'DisplayName', 'TGF-a');
xlabel('Time (hours)', 'FontSize', 11);
ylabel('Normalized Concentration (0-1)', 'FontSize', 11);
ylim([0, 1.1]);
legend('Location', 'best', 'FontSize', 10);
title('HIF-1a vs TGF-a (Normalized)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
hold off;

sgtitle('TGF-a Signaling Dynamics', 'FontSize', 14, 'FontWeight', 'bold');
%% =========================================================================
% FIGURE 2: Spatial Analysis (Last Timestep)
% =========================================================================

fprintf('Creating Figure 2: Spatial Analysis...\n');

% Get the final timestep (now properly sorted)
data_final = all_data{end};

fprintf('  Using final timestep: Time = %.1f hours, Cells = %d\n\n', data_final.Time(1), height(data_final));

fig2 = figure('Name', 'Spatial Analysis', 'Position', [100, 100, 1200, 800]);

% Plot 1: Oxygen gradient
subplot(2, 3, 1);
scatter(data_final.DistanceFromCenter, data_final.Oxygen, 50, data_final.CellState, 'filled', 'MarkerFaceAlpha', 0.6);
colorbar;
xlabel('Distance from Center (cell units)', 'FontSize', 11);
ylabel('Oxygen Concentration', 'FontSize', 11);
title('Oxygen Gradient', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 2: HIF-1a expression
subplot(2, 3, 2);
scatter(data_final.DistanceFromCenter, data_final.HIF1Alpha, 50, data_final.Oxygen, 'filled', 'MarkerFaceAlpha', 0.6);
colorbar;
xlabel('Distance from Center (cell units)', 'FontSize', 11);
ylabel('HIF-1a Expression', 'FontSize', 11);
title('HIF-1a vs Distance', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 3: TGF-a expression
subplot(2, 3, 3);
scatter(data_final.DistanceFromCenter, data_final.TGFAlpha, 50, data_final.Oxygen, 'filled', 'MarkerFaceAlpha', 0.6);
colorbar;
xlabel('Distance from Center (cell units)', 'FontSize', 11);
ylabel('TGF-a Expression', 'FontSize', 11);
title('TGF-a vs Distance', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 4: HIF-1a vs Oxygen
subplot(2, 3, 4);
scatter(data_final.Oxygen, data_final.HIF1Alpha, 60, data_final.DistanceFromCenter, 'filled', 'MarkerFaceAlpha', 0.6);
colorbar;
xlabel('Oxygen Level', 'FontSize', 11);
ylabel('HIF-1a Expression', 'FontSize', 11);
title('HIF-1a Activation', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 5: EGFR-TGF-a crosstalk
subplot(2, 3, 5);
scatter(data_final.TGFAlpha, data_final.EGFRActivation, 60, data_final.Oxygen, 'filled', 'MarkerFaceAlpha', 0.6);
colorbar;
xlabel('TGF-a Level', 'FontSize', 11);
ylabel('EGFR Activation', 'FontSize', 11);
title('EGFR-TGF-a Signaling', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 6: Cell density radial profile
subplot(2, 3, 6);
radii = 0:0.5:ceil(max(data_final.DistanceFromCenter));
density = [];
for r = radii
    count = sum(abs(data_final.DistanceFromCenter - r) < 0.25);
    density = [density, count];
end
bar(radii, density, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
xlabel('Distance from Center (cell units)', 'FontSize', 11);
ylabel('Cell Count', 'FontSize', 11);
title('Cell Spatial Distribution', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

sgtitle(sprintf('Spatial Analysis at Time = %.1f hours (N=%d cells)', data_final.Time(1), height(data_final)), ...
        'FontSize', 14, 'FontWeight', 'bold');

%% =========================================================================
% FIGURE 3: 2D Cell Position Visualization (Oxygen)
% =========================================================================

fprintf('Creating Figure 3: 2D Visualizations...\n');

fig3 = figure('Name', '2D Cell Distribution', 'Position', [100, 100, 1400, 900]);

% Select 4 timepoints
target_times = [0, 49, 87, 124];
timesteps_to_plot = [];
times_extracted = cellfun(@(x) x.Time(1), all_data);

for target_time = target_times
    [~, idx] = min(abs(times_extracted - target_time));
    timesteps_to_plot = [timesteps_to_plot, idx];
end

theta = linspace(0, 2*pi, 100);

for idx = 1:4
    subplot(2, 2, idx);
    
    t_idx = timesteps_to_plot(idx);
    data_plot = all_data{t_idx};
    
    % Plot cells colored by oxygen
    scatter(data_plot.PositionX, data_plot.PositionY, 50, data_plot.Oxygen, 'filled', 'MarkerFaceAlpha', 0.8);
    cbar = colorbar;
    set(cbar, 'Ticks', linspace(0.01, 0.20, 6));
    cbar.Label.String = 'Oxygen Level';
    
    % Draw boundary circle
    max_radius = max(data_plot.DistanceFromCenter) * 1.1;
    boundary_x = max_radius * cos(theta);
    boundary_y = max_radius * sin(theta);
    hold on;
    plot(boundary_x, boundary_y, 'k--', 'LineWidth', 2);
    
    axis equal;
    grid on;
    xlabel('X (cell units)', 'FontSize', 10);
    ylabel('Y (cell units)', 'FontSize', 10);
    caxis([0.01, 0.20]);
    title(sprintf('Time = %.1f h | %d cells', data_plot.Time(1), height(data_plot)), ...
          'FontSize', 11);
    colormap(gca, 'hot');
end

sgtitle('Spatial Cell Distribution - Oxygen Levels', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% =========================================================================
% FIGURE 4: Cell States in Space
% =========================================================================

fprintf('Creating Figure 4: Cell States...\n');

fig4 = figure('Name', 'Cell States', 'Position', [100, 100, 1400, 900]);

% Custom colormap: Blue (viable) -> Green (hypoxic) -> Red (necrotic)
cmap_states = [0, 0, 1;      % Blue (viable)
               0, 1, 0;      % Green (hypoxic)
               1, 0, 0];     % Red (necrotic)

for idx = 1:4
    subplot(2, 2, idx);
    
    t_idx = timesteps_to_plot(idx);
    data_plot = all_data{t_idx};
    
    % Plot cells colored by state
    scatter(data_plot.PositionX, data_plot.PositionY, 50, data_plot.CellState, 'filled', 'MarkerFaceAlpha', 0.8);
    colormap(gca, cmap_states);
    
    cbar = colorbar;
    set(cbar, 'Ticks', [1, 2, 3]);
    set(cbar, 'TickLabels', {'Viable', 'Hypoxic', 'Necrotic'});
    
    % Draw boundary
    max_radius = max(data_plot.DistanceFromCenter) * 1.1;
    boundary_x = max_radius * cos(theta);
    boundary_y = max_radius * sin(theta);
    hold on;
    plot(boundary_x, boundary_y, 'k--', 'LineWidth', 2);
    
    axis equal;
    grid on;
    xlabel('X (cell units)', 'FontSize', 10);
    ylabel('Y (cell units)', 'FontSize', 10);
    caxis([0.5, 3.5]);
    title(sprintf('Time = %.1f h | %d cells', data_plot.Time(1), height(data_plot)), ...
          'FontSize', 11);
end

sgtitle('Cell States: Viable (Blue) → Hypoxic (Green) → Necrotic (Red)', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% =========================================================================
% STEP 3: Save Results
% =========================================================================

fprintf('\nSaving analysis results...\n');

% Create output directory
output_dir = 'matlab_analysis_results/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Save figures
savefig(fig1, [output_dir 'Figure1_PopulationDynamics.fig']);
exportgraphics(fig1, [output_dir 'Figure1_PopulationDynamics.png'], 'Resolution', 150);

savefig(fig1b, [output_dir 'Figure1B_TGFa_Analysis.fig']);
exportgraphics(fig1b, [output_dir 'Figure1B_TGFa_Analysis.png'], 'Resolution', 150);

savefig(fig2, [output_dir 'Figure2_SpatialAnalysis.fig']);
exportgraphics(fig2, [output_dir 'Figure2_SpatialAnalysis.png'], 'Resolution', 150);

savefig(fig3, [output_dir 'Figure3_2DCellDistribution.fig']);
exportgraphics(fig3, [output_dir 'Figure3_2DCellDistribution.png'], 'Resolution', 150);

savefig(fig4, [output_dir 'Figure4_CellStates.fig']);
exportgraphics(fig4, [output_dir 'Figure4_CellStates.png'], 'Resolution', 150);

% Save time evolution data to CSV
T = table(time_data.time', time_data.total_cells', time_data.viable', ...
          time_data.hypoxic', time_data.necrotic', time_data.mean_oxygen', ...
          time_data.mean_hif1a', time_data.mean_tgfa', time_data.mean_egfr', ...
          time_data.spheroid_radius', ...
          'VariableNames', {'Time_hours', 'TotalCells', 'ViableCells', 'HypoxicCells', ...
                           'NecrotiCells', 'MeanOxygen', 'MeanHIF1a', 'MeanTGFa', ...
                           'MeanEGFR', 'SpheroidRadius'});

writetable(T, [output_dir 'population_dynamics.csv']);

% Save statistics
fid = fopen([output_dir 'statistics.txt'], 'w');
fprintf(fid, '===========================================\n');
fprintf(fid, 'PC9 Spheroid Simulation Statistics\n');
fprintf(fid, '===========================================\n\n');

fprintf(fid, 'SIMULATION PARAMETERS:\n');
fprintf(fid, '  Duration: %.1f hours (%.1f days)\n', time_data.time(end), time_data.time(end)/24);
fprintf(fid, '  Number of timepoints: %d\n\n', length(time_data.time));

fprintf(fid, 'INITIAL STATE (Time = %.1f hours):\n', time_data.time(1));
fprintf(fid, '  Total cells: %d\n', time_data.total_cells(1));

fprintf(fid, '\nFINAL STATE (Time = %.1f hours):\n', time_data.time(end));
fprintf(fid, '  Total cells: %d\n', time_data.total_cells(end));
fprintf(fid, '  Viable cells: %d (%.1f%%)\n', time_data.viable(end), ...
        100*time_data.viable(end)/time_data.total_cells(end));
fprintf(fid, '  Hypoxic cells: %d (%.1f%%)\n', time_data.hypoxic(end), ...
        100*time_data.hypoxic(end)/time_data.total_cells(end));
fprintf(fid, '  Necrotic cells: %d (%.1f%%)\n', time_data.necrotic(end), ...
        100*time_data.necrotic(end)/time_data.total_cells(end));
fprintf(fid, '  Spheroid radius: %.2f cell units\n\n', time_data.spheroid_radius(end));

fprintf(fid, 'MEAN CONCENTRATIONS (Final):\n');
fprintf(fid, '  Oxygen: %.4f\n', time_data.mean_oxygen(end));
fprintf(fid, '  HIF-1a: %.4f\n', time_data.mean_hif1a(end));
fprintf(fid, '  TGF-a: %.4f\n', time_data.mean_tgfa(end));
fprintf(fid, '  EGFR: %.4f\n\n', time_data.mean_egfr(end));

fprintf(fid, 'GROWTH STATISTICS:\n');
fprintf(fid, '  Population change: %.1f%%\n', ...
        100*(time_data.total_cells(end)-time_data.total_cells(1))/time_data.total_cells(1));
fprintf(fid, '  Mean growth rate: %.2f cells/hour\n', mean(growth_rate));
fprintf(fid, '  Max growth rate: %.2f cells/hour\n', max(growth_rate));

fprintf(fid, '\n===========================================\n');
fclose(fid);

fprintf('\n=== Analysis Complete ===\n');
fprintf('Results saved to: %s\n\n', output_dir);
fprintf('Files created:\n');
fprintf('  - Figure1_PopulationDynamics.fig/.png\n');
fprintf('  - Figure1B_TGFa_Analysis.fig/.png\n');
fprintf('  - Figure2_SpatialAnalysis.fig/.png\n');
fprintf('  - Figure3_2DCellDistribution.fig/.png\n');
fprintf('  - Figure4_CellStates.fig/.png\n');
fprintf('  - population_dynamics.csv\n');
fprintf('  - statistics.txt\n');
