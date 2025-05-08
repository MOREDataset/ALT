%% ** BoxPlot of ALT values distribution based on location and decade
clear
clc
file = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_preprocessedDataset.csv';
opts = detectImportOptions(file);
data = readtable(file, opts);

% Calculate the decade for each year
data.Decade = floor(data.Year/10)*10;

% Get unique decades and regions
decades = unique(data.Decade);
regions = unique(data.MainRegion);

% Create a new figure for the box chart
figure('units','normalized','OuterPosition',[0 0 1 1]);
hold on;

% Convert MainRegion and Decade to categorical
data.MainRegion = categorical(data.MainRegion);
data.Decade = categorical(data.Decade);

% Create the box chart with different colors for each decade
b = boxchart(data.MainRegion, data.ALT_Max, 'GroupByColor', data.Decade,'JitterOutliers','on','MarkerStyle','x');

% Add labels and title
% xlabel('Region','FontSize',20,'Interpreter','latex');
ylabel('ALT (cm)','FontSize',20,'Interpreter','latex');
% title('Distribution of ALT during Each Decade Grouped by Main Region');
legend(categories(data.Decade), 'Location', 'northoutside','Interpreter','latex','Orientation','horizontal');

hold off;
grid on;
set(gca, 'Box', 'on', 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% exportgraphics(gcf, 'ALT_boxplot_decades.pdf','Resolution',300)


%% ** Feature Distribution Per Main Region
clear;
clc;

% Load dataset
horizon = 0;
train_data = readtable(strjoin(["TrainDatawERA5_", num2str(horizon), "Years_V3.csv"],''));
test_data = readtable(strjoin(["TestDatawERA5_", num2str(horizon), "Years_V3.csv"],''));
combined_data = [train_data; test_data];
combined_data.MainRegion = categorical(combined_data.MainRegion);
combined_data.NDVI = combined_data.NDVI * 0.0001;


% Define feature categories and main regions
features = {'PFG_days','SFG_days','NFG_days','NDVI','total_precipitation_sum','u_component_of_wind_10m'};%data.Properties.VariableNames(4:end-1);
main_regions = unique(combined_data.MainRegion); % Extract unique regions
colors = lines(length(main_regions)); % Different colors for each region

% Extract main region column (update this based on actual column name)
region_labels = categorical(combined_data.MainRegion); % Convert to categorical variable

% Compute category means (modify based on actual function usage)
[~, ~, ~, ~, ~, Vegetation, PFG, SFG, NFG, ERA5] = computeCategoryMeans(combined_data,3);

% Create figure with violin plots
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('horizontal','TileSpacing','compact','Padding','tight','units','normalized','outerposition',[0 0 1 1]);

numCategories = length(features);

for i = 1:numCategories
    category = features{i};
    nexttile;
    hold on;
    
    % Extract category data
    category_data = combined_data.(category);
    
    % Convert to array and ensure proper indexing
    feature_values = table2array(category_data); % Convert table to matrix
    feature_values = feature_values(:); % Convert matrix to vector
    regions_numeric = repelem(region_labels, size(category_data, 2)); % Repeat labels for multiple features

    % Generate violin plot for feature distributions per region
    violinplot(feature_values, grp2idx(regions_numeric));

    % Customize plot
    title([category ' Features by Region'], 'FontSize', 20, 'Interpreter', 'latex');
    grid on;
    box on;
    set(gca, 'xticklabels', main_regions, 'YGrid', 'on', 'YMinorGrid', 'on', 'FontSize', 20, 'TickLabelInterpreter', 'latex');
    
    hold off;
end

xlabel(t, 'Main Region', 'FontSize', 20, 'Interpreter', 'latex');
ylabel(t, 'Feature Value Distribution', 'FontSize', 20, 'Interpreter', 'latex');

% Save figure
% exportgraphics(gcf, 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\FeatureDistributionByRegion.pdf', 'Resolution', 300);



%% **Violin Plot for Spearman Correlation Across Horizons**
clear;
clc;

% Define horizons and feature categories
horizons = [0, 1, 2, 5];
feature_categories = {'Vegetation', 'PFG', 'SFG', 'NFG', 'ERA5'};

% Initialize storage
results = struct();
p_values = struct();

for c = 1:length(feature_categories)
    results.(feature_categories{c}) = [];
    p_values.(feature_categories{c}) = [];
end

% Loop over horizons and compute correlations
for hIdx = 1:length(horizons)
    horizon = horizons(hIdx);
    
    % Load data
    train_data = readtable(sprintf('TrainDatawERA5_%dYears_V3.csv', horizon));
    test_data = readtable(sprintf('TestDatawERA5_%dYears_V3.csv', horizon));
    combined_data = [train_data; test_data];

    % Compute category means
    [~, ~, ~, ~, ~, Vegetation, PFG, SFG, NFG, ERA5] = computeCategoryMeans(combined_data,3);
    
    % Store correlations
    for cIdx = 1:length(feature_categories)
        category = feature_categories{cIdx};
        category_data = eval(category);
        
        % Compute Spearman correlation and p-value
        [corrs, pvals] = corr(table2array(category_data), combined_data.ALT_Max, ...
                              'Rows', 'complete', 'Type', 'Spearman');
        
        corrs_vec = corrs(:); % Flatten
        pvals_vec = pvals(:); % Flatten
        
        % Store results
        results.(category) = [results.(category); corrs_vec'];
        p_values.(category) = [p_values.(category); pvals_vec'];
    end
end


%% **Visualization with Violin Plots**
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('horizontal','TileSpacing','compact','Padding','tight');

numCategories = length(feature_categories);
for i = 1:numCategories
    category = feature_categories{i};
    nexttile;
    hold on;
    
    % Extract correlation data
    data_to_plot = results.(category);

    % Convert data into column format for violinplot
    ydata = data_to_plot(:);  % Flatten all data into a single vector
    xgroupdata = repelem(horizons, size(data_to_plot, 1))'; % Expand horizons for each value
    
    % Create violin plot
    violinplot(ydata, xgroupdata, 'ViolinAlpha', 0.6, 'ViolinColor', [0 0.447 0.741]);

    % Compute means per horizon
    means = mean(data_to_plot, 2);

    % Overlay mean values with markers
    plot(horizons, means, 'ro-', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
         'LineWidth', 2, 'DisplayName', 'Mean');

    % Customize subplot
    title([category ' Features'], 'FontSize', 16, 'Interpreter', 'latex');
    grid on;
    box on;
    ylim([-0.3, 0.3]);
    xticks(horizons);
    xticklabels({'+0', '+1', '+2', '+5'});
    set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');

    hold off;
end

% Shared axis labels
xlabel(t, 'Prediction Horizon (Years)', 'FontSize', 18, 'Interpreter', 'latex');
ylabel(t, 'Spearman Correlation Coefficient ($\rho$)', 'FontSize', 18, 'Interpreter', 'latex');

% Save plot
% exportgraphics(gcf, 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_CorrOverHorizons_cat.pdf', 'Resolution', 300);

%% ** Spearman Correlation plot (input categories vs ALT)
clear;
clc;

% Define horizons and feature categories
horizons = [0, 1, 2, 5];
feature_categories = {'Vegetation', 'PFG', 'SFG', 'NFG', 'ERA5'};
colors = lines(length(horizons)); % Different colors for scatter points

% Initialize storage for correlation results
results = struct();
for c = 1:length(feature_categories)
    results.(feature_categories{c}) = [];
end

% Loop over horizons and compute correlations
for hIdx = 1:length(horizons)
    horizon = horizons(hIdx);
    
    % Load data
    train_data = readtable(sprintf('TrainDatawERA5_%dYears_V3.csv', horizon));
    test_data = readtable(sprintf('TestDatawERA5_%dYears_V3.csv', horizon));
    combined_data = [train_data; test_data];
    
    % Compute category means (modify based on actual function usage)
    [~, ~, ~, ~, ~, Vegetation, PFG, SFG, NFG, ERA5] = computeCategoryMeans(combined_data,3);
    
    % Store correlations
    for cIdx = 1:length(feature_categories)
        category = feature_categories{cIdx};
        category_data = eval(category);
        
        % Compute Spearman correlation
        corrs = corr(table2array(category_data), combined_data.ALT_Max, ...
                     'Rows', 'complete', 'Type', 'Spearman');
        corrs_vec = corrs(:); % Flatten

        % Store correlations
        results.(category) = [results.(category); corrs_vec'];
    end
end


% Create subplots for each category
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('horizontal','TileSpacing','compact','Padding','tight','units','normalized','outerposition',[0 0 1 1]);

numCategories = length(feature_categories);

for i = 1:numCategories
    category = feature_categories{i};
    nexttile,
    hold on;
    
    % Extract correlation data
    data_to_plot = results.(category);
    
    % Compute min and max for each horizon
    min_vals = min(data_to_plot, [], 2);
    max_vals = max(data_to_plot, [], 2);
    means = mean(data_to_plot, 2);
    
    % Compute error bar values (upper and lower range)
    lower_error = means - min_vals;
    upper_error = max_vals - means;
    
    % Plot error bars (showing min and max)
    er = errorbar(categorical(horizons), means, lower_error, upper_error, 'vertical', ...
        'LineStyle', 'none', 'Color', 'k', 'CapSize', 10,'LineWidth', 1.5, 'HandleVisibility','off');
    plot(categorical(horizons), means,'LineStyle', ':','Marker', 'o', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerEdgeColor', 'k','LineWidth', 1.5, 'DisplayName', 'Average value')
    
    % Customize subplot
    title([category ' features'], 'FontSize', 20, 'Interpreter','latex');
    grid on;
    box on;
    ylim([-.3 .31]);
    % if ~isequal(i,1)
    %     set(gca,'yticklabels',[])
    % end
    set(gca, 'xticklabels', {'+0', '+1', '+2', '+5'},'YGrid','on','YMinorGrid','on', 'FontSize', 20, 'TickLabelInterpreter', 'latex');
    hold off;
end
xlabel(t, 'Prediction horizon (years)', 'FontSize', 20, 'interpreter','latex');
ylabel(t, 'Spearman correlation coefficient ($\rho$)', 'FontSize', 20, 'interpreter','latex');
% Set shared y-axis properties
for i=1:4
    set(t.Children(i),'yticklabels',[])
end

lgd = legend('show', 'Location', 'northoutside');
lgd.Layout.Tile = 'north'; % Place the legend above the plots
set(lgd,'Interpreter','latex','FontSize',20);
% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_CorrOverHorizons_cat.pdf','Resolution',300)


%% ** Line Plot: ALT Sequences Over Time
clear
clc
% Load the dataset
filename = "ALT_preprocessedDataset.csv";
data = readtable(filename);

% Convert Year and MainRegion to categorical
data.Year = categorical(data.Year);
data.MainRegion = categorical(data.MainRegion);

% Unique regions and their styles
regions = categories(data.MainRegion);
markers = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'}; % Cycle markers if needed
regionColors = parula(numel(regions));
regionColors(end-1:end,:) = [0.9290 0.6940 0.1250;0.6350 0.0780 0.1840];
regionColors(3,:) = [0.2422    0.1504    0.6603];
regionColors(1,:) = [0.1085    0.6669    0.8734];
numMarkers = numel(markers);

% Track legend status per region
legendShown = containers.Map(regions, false(size(regions)));

% Plot setup
figure('units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(1, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
ax = nexttile;
hold(ax, 'on');

% Unique locations
locations = unique([data.Latitude, data.Longitude], 'rows');

% Plot all locations with region-specific style
for i = 1:size(locations, 1)
    loc_data = data(data.Latitude == locations(i,1) & data.Longitude == locations(i,2), :);
    region = string(loc_data.MainRegion(1));
    regionIdx = find(strcmp(regions, region));

    % Only display legend once per region
    if ~legendShown(region)
        displayName = region;
        legendShown(region) = true;
        vis = 'on';
    else
        displayName = ''; % Skip legend entry
        vis = 'off';
    end

    plot(ax, loc_data.Year, loc_data.ALT_Max, '-', ...
        'Color', [regionColors(regionIdx, :) 0.35], ...
        'Marker', markers{mod(regionIdx-1, numMarkers)+1}, ...
        'MarkerFaceColor', regionColors(regionIdx, :), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerSize', 6, ...
        'LineWidth', 1.5, ...
        'DisplayName', displayName, 'HandleVisibility',vis);
end

% Aesthetics
xlabel('', 'Interpreter','latex');
ylabel('Maximum Active Layer Thickness (cm)', 'Interpreter','latex', 'FontSize', 30);
grid on;
set(ax, 'XTickLabelRotation', 90, ...
    'Box', 'on', ...
    'FontSize', 25, ...
    'TickLabelInterpreter', 'latex', ...
    'YDir', 'reverse'); % Reverse y-axis

legend(ax, 'Location', 'northoutside', 'Interpreter','latex', 'FontSize',30, 'Orientation','horizontal');

set(gca, 'YDir', 'reverse');  % Reverse the current y-axis (which is the original x-axis)

exportgraphics(gca,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALTtrends_overtimeV2.pdf', 'Resolution', 300) % Use vector for best quality


%% ** Map CALM locations - Train/Test splits
close('all')
clear
clc
%load data
hrz=0;
test_file = sprintf('ALT_hrz%d_TestPredictions.csv', hrz);
train_file = sprintf('ALT_hrz%d_TrainPredictions.csv', hrz);
test_data = readtable(test_file);
train_data = readtable(train_file);

% figure('Units','normalized','Position',[0.05 0.2 0.9 0.5]);
figure('Position',  [203         353        1533         551])
gx = geoaxes("Basemap","grayland");
hold on;


% Plot Train sites
geoscatter(train_data.Latitude, train_data.Longitude, 50, ...
    [40 150 206]./255, 'v', 'filled', 'DisplayName', 'Trainset sites');

% Plot Temporal Test-only sites
geoscatter(test_data.Latitude, test_data.Longitude, 50, ...
    [188 104 124]./255, '^', 'filled', 'DisplayName', 'Testset sites');
% Adjust map view
geolimits([50 83], [-180 180]);
set(gca, 'FontSize', 28);

% Add legend
legend1 = legend(gca,'show');
set(legend1,'Interpreter','latex','FontSize',30,'Orientation','horizontal','Location','northoutside');
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strcat(string(gx.LatitudeAxis.TickValues), '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strcat(string(gx.LongitudeAxis.TickValues), '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
gx.LatitudeAxis.TickLabels = strcat(string(gx.LatitudeAxis.TickValues), '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strcat(string(gx.LongitudeAxis.TickValues), '$^{\circ}$');

exportgraphics(gca,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_DatasetDistributions_V2.pdf', 'Resolution', 300) % Use vector for best quality


%% **Plot of Chittarjee's Correlation of geospatial groups Across Horizons**
clear;
clc;

% Define horizons and feature categories
horizons = [0, 1, 2, 5]; % Numeric x-axis values
feature_categories = {'Vegetation', 'ERA5', 'PFG', 'SFG', 'NFG'};

% Initialize storage
results = struct();

for c = 1:length(feature_categories)
    results.(feature_categories{c}) = [];
end

% Loop over horizons and compute correlations
for hIdx = 1:length(horizons)
    horizon = horizons(hIdx);
    
    % Load and preprocess data
    train_data = readtable(sprintf('TrainDatawERA5_%dYears_V3.csv', horizon));
    test_data = readtable(sprintf('TestDatawERA5_%dYears_V3.csv', horizon));
    combined_data = [train_data; test_data];
    combined_data(:,4:end-1) = fillmissing(combined_data(:,4:end-1),"constant",0);
    combined_data.NDVI = combined_data.NDVI * 0.0001;
    combined_data.EVI = combined_data.EVI * 0.0001;

    % Compute category means
    [~, Coordinates, year, LC, Topography, Vegetation, PFG, SFG, NFG, ERA5] = computeCategoryMeans(combined_data, 3);
    
    % Store correlations
    for cIdx = 1:length(feature_categories)
        category = feature_categories{cIdx};
        category_data = eval(category);
        
        % Compute Chatterjee correlation
        corrs = chatterjee_corr_vectorized(category_data, combined_data.ALT_Max);
        corrs_vec = corrs(:); % Flatten

        % Append as new column (features × horizons)
        results.(category) = [results.(category), corrs_vec];
    end
end

% Plotting
figure('units','normalized','outerposition',[-0.0047    0.0833    0.8151    0.5167])
t = tiledlayout(1, length(feature_categories), 'TileSpacing', 'compact', 'Padding', 'tight');

for i = 1:length(feature_categories)
    category = feature_categories{i};
    nexttile;
    hold on;

    ydata = results.(category);  % (features × horizons)

    % Simulate swarm plot: scatter with horizontal jitter
    for h = 1:length(horizons)
        xvals = horizons(h) + 0.1 * randn(size(ydata,1),1); % jittered x
        scatter(xvals, ydata(:,h), 25, 'filled', 'MarkerFaceAlpha', 0.6);
    end

    % Plot horizon-wise means
    means = mean(ydata, 1);
    p = plot(horizons, means, 'ko-', 'MarkerSize', 8, ...
             'MarkerFaceColor', 'w', 'LineWidth', 1.5, ...
             'DisplayName', 'Average $\xi$');  % changed label to Chatterjee

    % Customize
    title([category ' Features'], 'FontSize', 18, 'Interpreter', 'latex');
    grid on;
    box on;
    ylim([0, 0.4]);
    xlim([-0.5, 5.5]);
    xticks(horizons);
    xticklabels({'+0', '+1', '+2', '+5'});
    set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');

    hold off;
end

% Shared labels
xlabel(t, 'Prediction Horizon (Years)', 'FontSize', 20, 'Interpreter', 'latex');
ylabel(t, '$\xi_n$', 'FontSize', 20, 'Interpreter', 'latex');

% % Legend (only once)
% lgd = legend(p(1), 'Location', 'northoutside');
% lgd.Layout.Tile = 'north';
% set(lgd, 'Interpreter', 'latex', 'FontSize', 16);

% Set shared y-axis properties
for i=1:4
    set(t.Children(i), 'yticklabels', []);
end
% Save plot
exportgraphics(gcf, 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_CorrOverHorizons_Violin.pdf', 'Resolution', 300);
%% ** Swarmchart Plot of Selected Features by Main Region
close('all');
clear;
clc;

% Load data for horizon 0
horizon = 0;
train_data = readtable(strjoin(["TrainDatawERA5_", num2str(horizon), "Years_V3.csv"],''));
test_data = readtable(strjoin(["TestDatawERA5_", num2str(horizon), "Years_V3.csv"],''));
data = [train_data; test_data];

% Ensure main region is categorical
data.MainRegion = categorical(data.MainRegion);
Y = data.MainRegion;
Ycat = unique(Y);

% Scale NDVI
data.NDVI = data.NDVI * 0.0001;

% Define selected features
features = {'PFG_days', 'SFG_days', 'NFG_days', 'skin_temperature', 'NDVI','Elevation'};
featurenames = {'PFG days', 'SFG days', 'NFG days', 'Skin Temperature', 'NDVI', 'Elevation'};
% Extract unique main regions
regions = unique(data.MainRegion);

% Create figure with flow layout
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'tight', 'units', 'normalized', 'outerposition', [0 0 1 1]);

for j = 1:length(features)
    nexttile;
    hold on;

    feature = features{j};
    % Extract feature data
    for ii = 1:length(Ycat)
        feature_values = data.(feature)(Y==Ycat(ii));
        if isequal(feature, 'skin_temperature')% remove missing
            feature_values(feature_values==0)=[];
        end
        % Generate violin plot
        swarmchart(categorical(strings(size(feature_values'))+string(Ycat(ii))), feature_values',10);
    end
    % Customize plot
    % xlabel('Main Region', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(strrep(featurenames{j}, '_', '\_'), 'Interpreter', 'latex', 'FontSize', 20);
    % title([strrep(feature, '_', ' '), ' Distribution'], 'FontSize', 20, 'Interpreter', 'latex');
    set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex', 'XGrid','on', 'XTickLabelRotation',90);
    if j < 4 
        set(gca, 'xticklabels', [])
    end
    box on;
    hold off;
    
end

% Save figure
exportgraphics(gcf, 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_SwarmFeatureByRegion.pdf', 'Resolution', 300);

%% **New figure: average ALT for each ecotype over the years
clear;
clc;

% Define horizons and feature categories
horizons = [0, 1, 2, 5]; % Numeric x-axis values
% feature_categories = ;

% Initialize storage
results = struct();

for c = 1:length(feature_categories)
    results.(feature_categories{c}) = [];
end

% Loop over horizons and compute correlations
for hIdx = 1:length(horizons)
    horizon = horizons(hIdx);
    
    % Load and preprocess data
    train_data = readtable(sprintf('TrainDatawERA5_%dYears_V3.csv', horizon));
    test_data = readtable(sprintf('TestDatawERA5_%dYears_V3.csv', horizon));
    combined_data = [train_data; test_data];
    combined_data(:,4:end-1) = fillmissing(combined_data(:,4:end-1),"constant",0);
    combined_data.NDVI = combined_data.NDVI * 0.0001;
    combined_data.EVI = combined_data.EVI * 0.0001;

    % specifically for LC1
    combined_data.LC_Type1



end

% Plotting
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1, length(feature_categories), 'TileSpacing', 'compact', 'Padding', 'tight');

for i = 1:length(feature_categories)
    category = feature_categories{i};
    nexttile;
    hold on;

    ydata = results.(category);  % (features × horizons)

    % Simulate swarm plot: scatter with horizontal jitter
    for h = 1:length(horizons)
        xvals = horizons(h) + 0.1 * randn(size(ydata,1),1); % jittered x
        scatter(xvals, ydata(:,h), 25, 'filled', 'MarkerFaceAlpha', 0.6);
    end

    % Plot horizon-wise means
    means = mean(ydata, 1);
    p = plot(horizons, means, 'ko-', 'MarkerSize', 8, ...
             'MarkerFaceColor', 'w', 'LineWidth', 1.5, ...
             'DisplayName', 'Average $\xi$');  % changed label to Chatterjee

    % Customize
    title([category ' Features'], 'FontSize', 16, 'Interpreter', 'latex');
    grid on;
    box on;
    ylim([0, 0.4]);
    xlim([-0.5, 5.5]);
    xticks(horizons);
    xticklabels({'+0', '+1', '+2', '+5'});
    set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');

    hold off;
end

% Shared labels
xlabel(t, 'Prediction Horizon (Years)', 'FontSize', 18, 'Interpreter', 'latex');
ylabel(t, 'Chatterjee Correlation Coefficient ($\xi$)', 'FontSize', 18, 'Interpreter', 'latex');

% % Legend (only once)
% lgd = legend(p(1), 'Location', 'northoutside');
% lgd.Layout.Tile = 'north';
% set(lgd, 'Interpreter', 'latex', 'FontSize', 16);

% Set shared y-axis properties
for i=1:4
    set(t.Children(i), 'yticklabels', []);
end
%% ** Naive Baselines and subcomponents of proposed model (as Table)
%%Compute performance of Naive Baselines (1-3)
clear
clc
load('ALT_TrainData_targets.mat');
load('ALT_TemporalTestData_targets.mat');
testData = temporalTestData;

%Define Prediction Horizons
horizons = [0, 1, 2, 5];

%Define Unique Locations
locations = unique(trainData(:, {'Lat', 'Long'}), 'rows');

results = {};

n_boot = 1000; % Number of bootstrap samples
rng(42); % For reproducibility

for h = horizons
    fprintf('📂 Processing Baseliens at horizon: %d\n', h)

    %Remove the latest `h` years from the training set
    maxYearPerLoc = groupsummary(trainData, {'Lat', 'Long'}, 'max', 'Year');
    trainDataH = trainData;

    for i = 1:height(maxYearPerLoc)
        lat_i = maxYearPerLoc.Lat(i);
        long_i = maxYearPerLoc.Long(i);
        maxYear = maxYearPerLoc.max_Year(i);
        trainDataH(trainDataH.Lat == lat_i & trainDataH.Long == long_i & ...
                   trainDataH.Year > (maxYear - h), :) = [];
    end

    %Location-Agnostic, Yearly Average Baseline
    yearlyMeanALT = groupsummary(trainDataH, 'Year', 'mean', 'ALT_Max');
    testData2 = outerjoin(testData, yearlyMeanALT, "Type", "left", 'Keys', 'Year', 'MergeKeys', true);
    predYearAvg = testData2.mean_ALT_Max;

    %Combined Region-Year Baseline
    regionYearMeanALT = groupsummary(trainDataH, {'MainRegion'}, 'mean', 'ALT_Max');
    testData2 = outerjoin(testData, regionYearMeanALT, "Type", "left", 'Keys', {'MainRegion'}, 'MergeKeys', true);
    predRegionAvg = testData2.mean_ALT_Max;

    %Time-Agnostic, Location Average Baseline
    locationMeanALT = groupsummary(trainDataH, {'Lat', 'Long'}, 'mean', 'ALT_Max');
    testData2 = outerjoin(testData, locationMeanALT, "Type", "left", 'Keys', {'Lat', 'Long'}, 'MergeKeys', true);
    predLocAvg = testData2.mean_ALT_Max;

    %Compute Metrics with CI
    actual = testData.ALT_Max;
    models = {'YearAvg','RegionAvg', 'LocAvg'};
    preds = {predYearAvg, predRegionAvg, predLocAvg};

    for i = 1:length(models)
        pred = preds{i};
        validIdx = ~isnan(pred) & ~isnan(actual);
        pred = pred(validIdx);
        actualValid = actual(validIdx);

        % Bootstrap estimates
        n = numel(actualValid);
        rmse_vals = zeros(n_boot,1);
        mae_vals = zeros(n_boot,1);
        r2_vals = zeros(n_boot,1);

        for b = 1:n_boot
            idx = randi(n, n, 1);
            y_true = actualValid(idx);
            y_pred = pred(idx);

            rmse_vals(b) = sqrt(mean((y_true - y_pred).^2));
            mae_vals(b) = mean(abs(y_true - y_pred));
            r2_vals(b) = 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2);
        end

        % Compute mean and CI
        RMSE = mean(rmse_vals);
        RMSE_CI = prctile(rmse_vals, [2.5 97.5]);

        MAE = mean(mae_vals);
        MAE_CI = prctile(mae_vals, [2.5 97.5]);

        R2 = mean(r2_vals);
        R2_CI = prctile(r2_vals, [2.5 97.5]);

        % Store Results
        results = [results; {
            h, models{i}, ...
            RMSE, RMSE_CI(1), RMSE_CI(2), ...
            MAE, MAE_CI(1), MAE_CI(2), ...
            R2, R2_CI(1), R2_CI(2)}];
    end
end

%Convert to Table and Display
results_baselines = cell2table(results, ...
    'VariableNames', {'Horizon', 'Method', ...
    'Test_RMSE', 'Test_RMSE_L', 'Test_RMSE_U', ...
    'Test_MAE', 'Test_MAE_L', 'Test_MAE_U', ...
    'Test_R2', 'Test_R2_L', 'Test_R2_U'});


%%Proposed model results
horizons = [0, 1, 2, 5];
methods = {'CatBoost', 'ExtraTrees', 'BaggingRegressor', 'ensemble'};
metrics = {'RMSE', 'MAE', 'R2'};

results = [];

for h = horizons
    fprintf('📂 Processing Proposed Model at Horizon: %d\n', h);
    test_file = sprintf('ALT_hrz%d_TestPredictions.csv', h);
    train_file = sprintf('ALT_hrz%d_TrainPredictions.csv', h);

    test_data = readtable(test_file);
    train_data = readtable(train_file);

    y_train = train_data.y_train;
    y_test = test_data.y_test;

    for m = 1:numel(methods)
        method = methods{m};

        % Train
        y_pred_train = train_data.(sprintf('y_trpred_%s', method));
        [rmse_tr, ci_rmse_tr] = compute_ci(y_train, y_pred_train, @rmse);
        [mae_tr, ci_mae_tr] = compute_ci(y_train, y_pred_train, @mae);
        [r2_tr, ci_r2_tr] = compute_ci(y_train, y_pred_train, @r2score);

        % Test
        y_pred_test = test_data.(sprintf('y_tepred_%s', method));
        [rmse_te, ci_rmse_te] = compute_ci(y_test, y_pred_test, @rmse);
        [mae_te, ci_mae_te] = compute_ci(y_test, y_pred_test, @mae);
        [r2_te, ci_r2_te] = compute_ci(y_test, y_pred_test, @r2score);

        % Store
        results = [results; {
            h, method, ...
            rmse_tr ci_rmse_tr(1), ci_rmse_tr(2), ...
            mae_tr, ci_mae_tr(1), ci_mae_tr(2), ...
            r2_tr, ci_r2_tr(1), ci_r2_tr(2), ...
            rmse_te, ci_rmse_te(1), ci_rmse_te(2), ...
            mae_te, ci_mae_te(1), ci_mae_te(2), ...
            r2_te, ci_r2_te(1), ci_r2_te(2)}];
    end
end

% Display Table
results_tbl = cell2table(results, ...
    'VariableNames', {'Horizon', 'Method', ...
    'Train_RMSE', 'Train_RMSE_L', 'Train_RMSE_U', ...
    'Train_MAE', 'Train_MAE_L', 'Train_MAE_U', ...
    'Train_R2', 'Train_R2_L', 'Train_R2_U', ...
    'Test_RMSE', 'Test_RMSE_L', 'Test_RMSE_U', ...
    'Test_MAE', 'Test_MAE_L', 'Test_MAE_U', ...
    'Test_R2', 'Test_R2_L', 'Test_R2_U'});

% results = [results_tbl(:,[1:2 12:end]); results_baselines]; % if the baselines are needed
results = results_tbl;
results = sortrows(results,"Horizon","ascend");

% Define the models to include in the report
models_to_display = {'CatBoost', 'ExtraTrees', 'BaggingRegressor', 'ensemble'};%'YearAvg', 'RegionAvg', 'LocAvg',
metric_names = {'Test_RMSE', 'Test_MAE', 'Test_R2'};% 'Train_RMSE', 'Train_R2',
metric_labels = {'RMSE', 'MAE', 'R^2'};%'MAE', 
horizons = unique(results.Horizon);

fprintf('\n📊 ALT Prediction Benchmark (with 95%% CI)\n\n');

% Print header
header = sprintf('%-30s', 'Model');
for h = horizons'
    for k = 1:numel(metric_labels)
        header = [header, sprintf('%8s (+%d)', metric_labels{k}, h)];
    end
end
fprintf('%s\n', header);
fprintf('%s\n', repmat('-', 1, length(header)));

% Loop over each model
for model = models_to_display
    line = sprintf('%s;', model{1});
    for h = horizons'
        row = strcmp(results.Method, model{1}) & results.Horizon == h;
        for m = 1:numel(metric_names)
            if any(row)
                mean_val = results{row, metric_names{m}};
                lb = results{row, [metric_names{m}, '_L']};
                ub = results{row, [metric_names{m}, '_U']};
                pm = (ub - lb) / 2;
                % line = [line, sprintf('%.2f±%.2f;', mean_val, pm)];
                line = [line, sprintf('%.3f;', mean_val)];
            else
                line = [line, sprintf('%s', '---')];
            end
        end
    end
    fprintf('%s\n', line);
end


%% ** Main Results Figure - part1
close('all')
clear
clc
%%Compute performance of Naive Baselines (1-3)
clear
clc
load('ALT_TrainData_targets.mat');
load('ALT_TemporalTestData_targets.mat');
testData = temporalTestData;

%Define Prediction Horizons
horizons = [0, 1, 2, 5];

%Define Unique Locations
locations = unique(trainData(:, {'Lat', 'Long'}), 'rows');

results = {};

n_boot = 1000; % Number of bootstrap samples
rng(42); % For reproducibility

for h = horizons
    fprintf('📂 Processing Baseliens at horizon: %d\n', h)

    %Remove the latest `h` years from the training set
    maxYearPerLoc = groupsummary(trainData, {'Lat', 'Long'}, 'max', 'Year');
    trainDataH = trainData;

    for i = 1:height(maxYearPerLoc)
        lat_i = maxYearPerLoc.Lat(i);
        long_i = maxYearPerLoc.Long(i);
        maxYear = maxYearPerLoc.max_Year(i);
        trainDataH(trainDataH.Lat == lat_i & trainDataH.Long == long_i & ...
                   trainDataH.Year > (maxYear - h), :) = [];
    end

    %Location-Agnostic, Yearly Average Baseline
    yearlyMeanALT = groupsummary(trainDataH, 'Year', 'mean', 'ALT_Max');
    testData2 = outerjoin(testData, yearlyMeanALT, "Type", "left", 'Keys', 'Year', 'MergeKeys', true);
    predYearAvg = testData2.mean_ALT_Max;
    % 
    % %Combined Region-Year Baseline
    % regionYearMeanALT = groupsummary(trainDataH, {'MainRegion'}, 'mean', 'ALT_Max');
    % testData2 = outerjoin(testData, regionYearMeanALT, "Type", "left", 'Keys', {'MainRegion'}, 'MergeKeys', true);
    % predRegionAvg = testData2.mean_ALT_Max;

    %Time-Agnostic, Location Average Baseline
    locationMeanALT = groupsummary(trainDataH, {'Lat', 'Long'}, 'mean', 'ALT_Max');
    testData2 = outerjoin(testData, locationMeanALT, "Type", "left", 'Keys', {'Lat', 'Long'}, 'MergeKeys', true);
    predLocAvg = testData2.mean_ALT_Max;

    % Permafrost Type-Based Baseline
    fprintf(' → Assigning permafrost type for test sites...\n')
    testCats = assignEnvironmentalCategories(testData.Lat, testData.Long);
    trainCats = assignEnvironmentalCategories(trainDataH.Lat, trainDataH.Long);
    testData.PermafrostType = testCats.permafrostType;
    trainDataH.PermafrostType = trainCats.permafrostType;
    
    % Average by Permafrost Type
    typeMeanALT = groupsummary(trainDataH, 'PermafrostType', 'mean', 'ALT_Max');
    testData2 = outerjoin(testData, typeMeanALT, "Type", "left", ...
        'Keys', 'PermafrostType', 'MergeKeys', true);
    predPermafrostAvg = testData2.mean_ALT_Max;

    %Compute Metrics with CI
    actual = testData.ALT_Max;
    models = {'YearAvg','LocAvg', 'PermafrostAvg'};
    preds = {predYearAvg, predLocAvg, predPermafrostAvg};

    for i = 1:length(models)
        pred = preds{i};
        validIdx = ~isnan(pred) & ~isnan(actual);
        pred = pred(validIdx);
        actualValid = actual(validIdx);

        % Bootstrap estimates
        n = numel(actualValid);
        rmse_vals = zeros(n_boot,1);
        mae_vals = zeros(n_boot,1);
        r2_vals = zeros(n_boot,1);

        for b = 1:n_boot
            idx = randi(n, n, 1);
            y_true = actualValid(idx);
            y_pred = pred(idx);

            rmse_vals(b) = sqrt(mean((y_true - y_pred).^2));
            mae_vals(b) = mean(abs(y_true - y_pred));
            r2_vals(b) = 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2);
        end

        % Compute mean and CI
        RMSE = mean(rmse_vals);
        RMSE_CI = prctile(rmse_vals, [2.5 97.5]);

        MAE = mean(mae_vals);
        MAE_CI = prctile(mae_vals, [2.5 97.5]);

        R2 = mean(r2_vals);
        R2_CI = prctile(r2_vals, [2.5 97.5]);

        % Store Results
        results = [results; {
            h, models{i}, ...
            RMSE, RMSE_CI(1), RMSE_CI(2), ...
            MAE, MAE_CI(1), MAE_CI(2), ...
            R2, R2_CI(1), R2_CI(2)}];
    end
end

%Convert to Table and Display
results_baselines = cell2table(results, ...
    'VariableNames', {'Horizon', 'Method', ...
    'Test_RMSE', 'Test_RMSE_L', 'Test_RMSE_U', ...
    'Test_MAE', 'Test_MAE_L', 'Test_MAE_U', ...
    'Test_R2', 'Test_R2_L', 'Test_R2_U'});


%%Plotting per horizon
% Define models and horizons
%%Proposed model results
methods = {'ensemble'};
metrics = {'RMSE', 'MAE', 'R2'};

results = [];

for h = [0, 1, 2, 5]
    fprintf('📂 Processing Proposed Model at Horizon: %d\n', h);
    test_file = sprintf('ALT_hrz%d_TestPredictions.csv', h);
    train_file = sprintf('ALT_hrz%d_TrainPredictions.csv', h);

    test_data = readtable(test_file);
    train_data = readtable(train_file);

    y_train = train_data.y_train;
    y_test = test_data.y_test;

    for m = 1:numel(methods)
        method = methods{m};

        % Train
        y_pred_train = train_data.(sprintf('y_trpred_%s', method));
        [rmse_tr, ci_rmse_tr] = compute_ci(y_train, y_pred_train, @rmse);
        [mae_tr, ci_mae_tr] = compute_ci(y_train, y_pred_train, @mae);
        [r2_tr, ci_r2_tr] = compute_ci(y_train, y_pred_train, @r2score);

        % Test
        y_pred_test = test_data.(sprintf('y_tepred_%s', method));
        [rmse_te, ci_rmse_te] = compute_ci(y_test, y_pred_test, @rmse);
        [mae_te, ci_mae_te] = compute_ci(y_test, y_pred_test, @mae);
        [r2_te, ci_r2_te] = compute_ci(y_test, y_pred_test, @r2score);

        % Store
        results = [results; {
            h, method, ...
            rmse_tr ci_rmse_tr(1), ci_rmse_tr(2), ...
            mae_tr, ci_mae_tr(1), ci_mae_tr(2), ...
            r2_tr, ci_r2_tr(1), ci_r2_tr(2), ...
            rmse_te, ci_rmse_te(1), ci_rmse_te(2), ...
            mae_te, ci_mae_te(1), ci_mae_te(2), ...
            r2_te, ci_r2_te(1), ci_r2_te(2)}];
    end
end

% Display Table
tbl = cell2table(results, ...
    'VariableNames', {'Horizon', 'Method', ...
    'Train_RMSE', 'Train_RMSE_L', 'Train_RMSE_U', ...
    'Train_MAE', 'Train_MAE_L', 'Train_MAE_U', ...
    'Train_R2', 'Train_R2_L', 'Train_R2_U', ...
    'Test_RMSE', 'Test_RMSE_L', 'Test_RMSE_U', ...
    'Test_MAE', 'Test_MAE_L', 'Test_MAE_U', ...
    'Test_R2', 'Test_R2_L', 'Test_R2_U'});

tbl = [tbl(:,[1:2 12:end]);results_baselines];
methods = unique(tbl.Method);
%% ** Main Results Figure - part2 
% Colors
colors = [0.2 0.6 1; 0.4660 0.6740 0.1880;0.9290 0.6940 0.1250;0.6350 0.0780 0.1840];
% colors = [
%     0.4660 0.6740 0.1880;
%     0.3010 0.7450 0.9330;
%     0.8500 0.3250 0.0980;
%     0.4940 0.1840 0.5560;
%     0.9290 0.6940 0.1250;
%     0.6350 0.0780 0.1840;
% ];

tbl.Method = categorical(tbl.Method, methods, 'Ordinal', true);
horizons = unique(tbl.Horizon);
hz_labels = compose('+%d days', horizons);
metrics = {'Test_RMSE', 'Test_MAE', 'Test_R2'};
metric_labels = {'RMSE (cm)', 'MAE (cm)', 'R$^2$'};

figure('Units','normalized','OuterPosition',[0 0 1 1]);
t = tiledlayout(3, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for k = 1:length(metrics)
    nexttile;
    hold on;

    data_vals = NaN(length(horizons), length(methods));
    data_err_low = NaN(length(horizons), length(methods));
    data_err_up = NaN(length(horizons), length(methods));

    for h = 1:length(horizons)
        for m = 1:length(methods)
            row = strcmp(string(tbl.Method), methods{m}) & tbl.Horizon == horizons(h);
            val = tbl{row, metrics{k}};
            lwr = tbl{row, [metrics{k}, '_L']};
            upr = tbl{row, [metrics{k}, '_U']};

            data_vals(h, m) = val;
            data_err_low(h, m) = lwr;
            data_err_up(h, m) = upr;
        end
    end

    b = bar(data_vals, 'grouped', 'FaceAlpha', 0.7);
    for m = 1:length(methods)
        b(m).FaceColor = colors(m, :);
        b(m).EdgeColor = 'none';
        b(m).DisplayName = methods{m};
    end
    legendHandles = b; % this array is in the correct order of bar series

    % Add asymmetric error bars
    [ngroups, nbars] = size(data_vals);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
        nmetrics = 3;
        nbars_per_metric = nbars / nmetrics;
        
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, data_vals(:,i), (data_err_up(:,i) - data_err_low(:,i))/2, ...
                'k', 'linestyle', 'none', 'LineWidth', 1.2, 'CapSize', 10);
           if i==4
               xregion(x-.09,x+.09,FaceColor=[0.941176470588235 0.941176470588235 0.941176470588235],EdgeColor='w')
           end
            % 🔍 Determine which metric this bar belongs to
            metric_idx = ceil(i / nbars_per_metric);
            metric_start = (metric_idx - 1) * nbars_per_metric + 1;
            metric_end = metric_start + nbars_per_metric - 1;
        
            % ✅ Compute max value **within this metric block only**
            curr_max = max(max(data_vals(:, metric_start:metric_end) + ...
                               data_err_up(:, metric_start:metric_end)));
        
            for j = 1:ngroups
                rot = 90;
                val = data_vals(j, i);
                bar_err = (data_err_up(j, i) - data_err_low(j, i)) / 2;
                ratio = val / curr_max;
        
                label = sprintf(['%.3f' char(10) '$\\pm$%.2f'], val, bar_err);
                x_shift = 0;
                ha = 'center';

                if i==4 && k~=3 %only when it is the proposed model%ratio < 0.46  && 
                    yloc = val + bar_err + 0.03 * curr_max;
                    color = 'k';
                    ha = 'left';
                    x_shift = 0;%0.03;
                    label = strrep(label, char(10), ''); % always multiline
                elseif ratio > 0.48 && k~=3
                    yloc = val / 2;
                    color = 'w';  % 👈 This now only applies in mid-range values
                    label = strrep(label, char(10), ''); % always multiline
                elseif k==2 && i==1 
                    yloc = val + bar_err + 0.03 * curr_max;
                    color = 'k';
                    ha = 'left';
                    x_shift = 0;%0.03;
                    label = strrep(label, char(10), ''); % always multiline
                elseif val<0 
                    yloc = 0.04;
                    color = 'k';
                    ha = 'left';
                    x_shift = 0;%0.03;                    
                    label = '$< 0$'; % don't report values if negative!
                    rot = 0;
                elseif [k~=3 && i~=4]
                    yloc = val / 2;
                    color = 'w';
                elseif k==3 && i==4
                    yloc = val / 2;
                    color = 'w';
                    ha = 'center';
                    x_shift = 0;%0.03;                    
                    label = strrep(label, char(10), ''); % always multiline
                elseif ratio>0.3
                    yloc = val /2;
                    color = 'w';
                    ha = 'center';
                    % label = strrep(label, char(10), ''); % always multiline
                else
                    yloc = val + bar_err + 0.03 * curr_max;
                    color = 'k';
                    ha = 'left';
                    label = strrep(label, char(10), ''); % always multiline
                end
        
                text(x(j) + x_shift, yloc, label, ...
                    'HorizontalAlignment', ha, ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 18, ...
                    'Interpreter', 'latex', ...
                    'Color', color, ...
                    'Rotation', rot);
            end
        end


    set(gca, 'XTick', 1:length(hz_labels), 'XTickLabel', [], ...
        'FontSize', 20, 'TickLabelInterpreter','latex');
    ylabel(metric_labels{k}, 'Interpreter', 'latex', 'FontSize', 20);
    grid on; box on;

    if k == 3
        ylim([0 1]);
        set(gca, 'XTick', 1:length(hz_labels), 'XTickLabel', hz_labels, ...
            'FontSize', 20, 'TickLabelInterpreter','latex');
    end
end

% lgd = legend({'Base Model (CatBoost)', 'Base Model + Correction', 'Proposed model'}, 'Location', 'northoutside', 'Orientation', 'horizontal');
% lgd = legend(legendHandles, methods, ...
%     'Location', 'northoutside', ...
%     'Orientation', 'horizontal', ...
%     'Interpreter', 'latex', ...
%     'FontSize', 20);
lgd = legend({'Naive Model-Location Average', 'Naive Model-Permafrost Average', 'Naive Model-Yearly Average', 'Proposed Model'}, 'Location', 'northoutside', 'Orientation', 'horizontal');
lgd.Layout.Tile = 'north';
set(lgd, 'Interpreter', 'latex', 'FontSize', 20);

exportgraphics(t,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_Barplot_BaselineBench.pdf'),'Resolution',300)

%% ** Results scatter plot per MainRegion for the Proposed Model (Alaska, Russia, other)
close all
clear
clc

% Target regions: Alaska, Russia, and "other" = everything else
target_regions = {'Alaska', 'Russia'};

figure('Units', 'normalized', 'OuterPosition', [0 0 0.7 1]);
t = tiledlayout(2, 3, 'TileSpacing', 'tight', 'Padding', 'tight');

% Parameters 
markers = {'o', 'd', 'v'};
colors = [0.2 0.6 1; 0.4660 0.6740 0.1880;0.6350 0.0780 0.1840];%0.9290 0.6940 0.1250;

for row = 1:2 % 1 for H=0, 2 for H=5
    h = (row - 1) * 5;

    % Load data
    test_file = sprintf('ALT_hrz%d_TestPredictions.csv', h);
    train_file = sprintf('ALT_hrz%d_TrainPredictions.csv', h);
    test_data = readtable(test_file);
    train_data = readtable(train_file);

    % Create a logical mask for the 'other' category
    is_other_test = ~ismember(test_data.MainRegion, target_regions);
    is_other_train = ~ismember(train_data.MainRegion, target_regions);

    % Regions to loop over (2 main + 1 combined "other")
    regions = [target_regions, {'other'}];
    region_masks_test = {
        strcmp(test_data.MainRegion, 'Alaska'), ...
        strcmp(test_data.MainRegion, 'Russia'), ...
        is_other_test
        };
    region_masks_train = {
        strcmp(train_data.MainRegion, 'Alaska'), ...
        strcmp(train_data.MainRegion, 'Russia'), ...
        is_other_train
        };

    for r = 1:3
        region = regions{r};
        idx_te = region_masks_test{r};
        idx_tr = region_masks_train{r};

        y_te = test_data.y_test(idx_te);
        y_te_pred = test_data.y_tepred_ensemble(idx_te);

        y_tr = train_data.y_train(idx_tr);
        y_tr_pred = train_data.y_trpred_ensemble(idx_tr);

        % Compute metrics
        [rmse_tr, ci_rmse_tr] = compute_ci(y_tr, y_tr_pred, @rmse);
        [mae_tr, ci_mae_tr] = compute_ci(y_tr, y_tr_pred, @mae);
        [r2_tr, ci_r2_tr] = compute_ci(y_tr, y_tr_pred, @r2score);
        [rmse_te, ci_rmse_te] = compute_ci(y_te, y_te_pred, @rmse);
        [mae_te, ci_mae_te] = compute_ci(y_te, y_te_pred, @mae);
        [r2_te, ci_r2_te] = compute_ci(y_te, y_te_pred, @r2score);

        % Scatter plot
        nexttile;
        hold on; grid on;

        s.(['p',num2str(r)]) = scatter(y_te, y_te_pred, 30, markers{r}, 'filled', ...
            'MarkerFaceColor', colors(r,:), 'MarkerFaceAlpha', 0.4, ...
            'MarkerEdgeColor', 'none');

        % Diagonal line
        minVal = min([y_te; y_te_pred]);
        maxVal = max([y_te; y_te_pred]);
        plot([minVal, maxVal], [minVal, maxVal], '--k', 'LineWidth', 1.2);

        set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');

        axis square; box on
        xlim([minVal, maxVal]);
        ylim([minVal, maxVal]);

        % Annotate metrics
        train_str = sprintf(['Training:' char(10) 'RMSE = %.2f $\\pm$ %.2f' char(10) ...
            'MAE = %.2f $\\pm$ %.2f' char(10) ...
            'R$^2$ = %.2f $\\pm$ %.2f'], ...
            rmse_tr, diff(ci_rmse_tr)/2, ...
            mae_tr, diff(ci_mae_tr)/2, ...
            r2_tr, diff(ci_r2_tr)/2);

        test_str = sprintf(['Testing:' char(10) 'RMSE = %.2f $\\pm$ %.2f' char(10) ...
            'MAE = %.2f $\\pm$ %.2f' char(10) ...
            'R$^2$ = %.2f $\\pm$ %.2f'], ...
            rmse_te, diff(ci_rmse_te)/2, ...
            mae_te, diff(ci_mae_te)/2, ...
            r2_te, diff(ci_r2_te)/2);

        text(maxVal - 0.05*(maxVal - minVal), minVal + 0.15*(maxVal - minVal), ...
            train_str, 'FontSize', 14, 'Interpreter', 'latex', ...
            'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

        text(minVal + 0.51*(maxVal-minVal), maxVal - 0.15*(maxVal-minVal), ...
            test_str, 'FontSize', 14, 'Interpreter', 'latex', ...
            'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

                % Title and formatting
        if [row == 1 && r == 1] || [row == 2 && r == 1]
            title(sprintf('Horizon = +%d year', h), 'BackgroundColor',[0.941176470588235 0.941176470588235 0.941176470588235], ...
                'FontSize', 18, 'Interpreter', 'latex', 'HorizontalAlignment','right', 'VerticalAlignment','bottom', Margin=.1);
        else
            set(gca, 'yticklabels', [])
        end
        if row ~= 2
            set(gca, 'xticklabels', [])
        end
        set(gca, 'XDir', 'reverse');  % Reverse the current y-axis (which is the original x-axis)
        set(gca, 'YDir', 'reverse');  % Reverse the current y-axis (which is the original x-axis)
    end
end

ylabel(t, 'Target ALT (cm)', 'FontSize', 20, 'Interpreter', 'latex');
xlabel(t, 'Predicted ALT (cm)', 'FontSize', 20, 'Interpreter', 'latex');
% Add a common legend
lgd = legend([s.p1,s.p2,s.p3], regions, 'Orientation', 'horizontal', 'Location', 'northoutside');
lgd.Layout.Tile = 'north'; % Place the legend above the plots
set(lgd,'Interpreter','latex','FontSize',20);

% export
exportgraphics(gcf, 'ALT_ScatterByRegion.pdf', 'Resolution', 300)

%% ** SpiderPlots - different parameters
addpath(genpath('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\SpiderPlot'))
clear;
clc;

categories = {'permafrostType', 'groundIceType', 'landformType'};
categories_text = {'Permafrost Type', 'Ground Ice Content', 'Landforms'};
category_labels_map = containers.Map;
category_labels_map('permafrostType') = {'Continuous', 'Discontinuous', 'Sporadic', 'Unknown'};
category_labels_map('groundIceType') = {'High', 'Low', 'Medium', 'Unknown'};
category_labels_map('landformType') = {'Lowlands', 'Mountains', 'Unknown'};

metrics = {'RMSE', 'MAE', 'R2'};
horizons = [0, 1, 2, 5];

% Load and tag all data once
allData = [];
for h = horizons
    data = readtable(sprintf('ALT_hrz%d_TestPredictions.csv', h));
    latitudes = data.Latitude;
    longitudes = data.Longitude;
    categoryTable = assignEnvironmentalCategories(latitudes, longitudes);
    categoryTable.permafrostType = categorical(categoryTable.permafrostType);
    categoryTable.groundIceType = categorical(categoryTable.groundIceType);
    categoryTable.landformType = categorical(categoryTable.landformType);
    data = [data categoryTable(:, 3:end)];
    data.horizon = repmat(h, height(data), 1);
    allData = [allData; data];
end
%
clc
fprintf('Category, horizon, RMSE, MAE, R2\n')
% Create radar plots
figure('units','normalized','OuterPosition',[0           0        .6         1]);
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact','units','normalized','outerposition',[0 0 1 1]);
for c = 1:length(categories)
    cat = categories{c};
    cat_labels = category_labels_map(cat);
    uniqueCats = string(cat_labels);
    
    % Preallocate matrix for metrics
    RMSE = nan(length(uniqueCats), length(horizons));
    MAE = nan(length(uniqueCats), length(horizons));
    R2 = nan(length(uniqueCats), length(horizons));
    
    for k = 1:length(uniqueCats)
        currentCategory = uniqueCats(k);
        categoryData = allData(allData.(cat) == currentCategory, :);

        for hIdx = 1:length(horizons)
            h = horizons(hIdx);
            thisH = categoryData(categoryData.horizon == h, :);
            y_true = thisH.y_test;
            y_pred = thisH.y_tepred_ensemble;

            if length(y_true) < 2
                RMSE(k, hIdx) = NaN;
                MAE(k, hIdx) = NaN;
                R2(k, hIdx) = NaN;
            else
                [rmse_val, ~] = compute_ci(y_true, y_pred, @rmse);
                [mae_val, ~]  = compute_ci(y_true, y_pred, @mae);
                [r2_val, ~]   = compute_ci(y_true, y_pred, @r2score);
                RMSE(k, hIdx) = rmse_val;
                MAE(k, hIdx) = mae_val;
                R2(k, hIdx) = r2_val;
                fprintf('- %s-%s, %d, %.2f, %.2f, %.2f\n', cat, currentCategory, h, rmse_val, mae_val, r2_val)
            end
        end
    end

    % Plot RMSE
    nexttile
    spider_plot(RMSE', ...
        'AxesLabels', cellstr(uniqueCats), ...
        'AxesInterval', 2, ...
        'AxesLimits', [zeros(1, height(RMSE)); ones(1, height(RMSE))*33],... % [min axes limits; max axes limits]
        'FillOption', {'on'}, ...
        'FillTransparency', 0.15, ...
        'LineWidth', 2, ...
        'AxesInterpreter', {'latex'},...
        'AxesTickInterpreter', {'latex'},...
        'AxesLabelsEdge', 'none','AxesPrecision', [0],...
        'AxesDisplay', 'one', 'AxesFontSize', 18, 'LabelFontSize', 18, 'BackgroundColor', 'none');
    % title(sprintf('%s', categories_text{c}), 'FontSize', 18, 'Interpreter','latex');
    % Plot R²
    nexttile
    % R2(isinf(R2)) = 0.5;
    spider_plot(R2', ...
        'AxesLabels', cellstr(uniqueCats), ...
        'AxesInterval', 2, ...
        'AxesLimits', [.5*ones(1, height(R2)); ones(1, height(R2))],... % [min axes limits; max axes limits]
        'FillOption', {'on'}, ...
        'FillTransparency', 0.15, ...
        'LineWidth', 2, ...
        'AxesInterpreter', {'latex'},...
        'AxesTickInterpreter', {'latex'},...
        'AxesLabelsEdge', 'none','AxesPrecision', [1],...
        'AxesDisplay', 'one', 'AxesFontSize', 18, 'LabelFontSize', 18, 'BackgroundColor', 'none');
end
legend1 = legend(gca,{'+0 year', '+1 year', '+2 years', '+5 years'}, Interpreter="latex", FontSize=18, Orientation="horizontal");
legend1.Layout.Tile = 'North'; % <-- place legend east of tiles
% exportgraphics(gcf,['C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_Results_spiderplot.pdf'], 'Resolution', 300) % Use vector for best quality    

%% ** Boxplots - Benchmark with other models
clear;
clc;

% Model names
models = {
    'EN', 'DT', 'RF', 'KNN', 'SVR', 'MLP', ...
    'ResNet50', 'DenseNet121', 'AutoInt', 'FTTransformer',...
    'Proposed Model'};

metric_names = {'RMSE (cm)', 'MAE (cm)', 'R$^2$'};

% Soft, friendly colors (color palette)
cmap = turbo(numel(models)); % Try also: parula, cool, autumn, etc.

% Metrics matrix (models x [RMSE0 MAE0 R²0 RMSE1 ...])
metrics = [
54.177	38.816	0.067	54.250	38.891	0.065	54.293	38.989	0.063	54.226	38.872	0.066
33.774	19.029	0.638	33.842	19.261	0.636	29.761	16.995	0.719	34.433	17.836	0.623
26.474	14.963	0.777	26.084	14.939	0.784	26.257	14.271	0.781	26.728	14.830	0.773
33.002	20.572	0.654	33.558	20.531	0.642	34.377	20.035	0.624	33.717	19.770	0.639
53.003	34.557	0.107	53.181	34.718	0.101	53.403	34.846	0.094	53.371	34.795	0.095
30.882	21.789	0.697	31.399	21.698	0.687	30.502	21.424	0.704	30.125	20.573	0.712
40.963	29.728	0.467	41.492	28.479	0.453	42.556	30.349	0.425	43.317	31.079	0.404
40.114	31.291	0.489	38.354	27.897	0.533	44.230	31.287	0.378	39.872	28.111	0.495
49.103	32.072	0.234	37.922	24.357	0.543	44.423	29.530	0.373	49.031	30.845	0.236
41.127	26.556	0.463	47.518	29.447	0.283	48.942	31.854	0.239	51.899	32.236	0.144
24.095	14.633	0.811	24.181	14.657	0.809	24.106	13.853	0.81	24.418	14.753	0.806
];

figure('Units','normalized','OuterPosition',[0 0 0.58 1]);
tiledlayout(3,1,'TileSpacing','tight','Padding','compact');

clrz = parula(length(models));

for m = 1:3
    nexttile;

    % Extract 4 values per model for the current metric
    data = zeros(numel(models)*4, 1);
    group = strings(numel(models)*4, 1);

    for i = 1:numel(models)
        for h = 1:4
            data((i-1)*4 + h) = metrics(i, (h-1)*3 + m);
            group((i-1)*4 + h) = models{i};
        end
    end

    % Plot
    boxplot(data, group, ...
        'Widths', 0.6, 'Whisker', 1.5);

    % Hide visible boxes
    set(findobj(gca, 'Tag', 'Box'), 'Visible', 'off');

    % Fancy colored boxes
    h = findobj(gca, 'Tag', 'Box');
    for j = 2:numel(h)
        patch(get(h(j),'XData'), get(h(j),'YData'), [41 139 195]./255, ...
            'FaceAlpha', 0.7, 'EdgeColor', 'none');
    end
    patch(get(h(1),'XData'), get(h(1),'YData'), [0.6350 0.0780 0.1840], ...
    'FaceAlpha', 0.7, 'EdgeColor', 'none');

    % Grid and styling
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 13;
    ylabel(metric_names{m}, 'Interpreter', 'latex', 'FontSize', 18);
    if m ~= 3
        set(gca, 'xticklabel', []);
        ylim([0 60])
        set(gca, 'YTick',[0:10:60])
    else
        ylim([0 1])
        set(gca, 'YTick',[0,0.2,0.4,0.6,0.8,1])
    end

    % Highlight median lines
    lines = findobj(gca, 'Tag', 'Median');
    set(lines(2:end), 'LineWidth', 1, 'Color', [0 0 0 1], 'LineStyle', '-');    
    xregion(10.7,11.3,FaceColor=[0.941176470588235 0.941176470588235 0.941176470588235],EdgeColor='w');
    set(lines(1), 'LineWidth', 1, 'Color', [0 0 0 1], 'LineStyle', '-');    

    set(gca, 'Box', 'on', 'FontSize', 20, 'TickLabelInterpreter', 'latex', ...
'YMinorGrid', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XTickLabelRotation',90);    
end
exportgraphics(gcf,['C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_ModelBenchmark.pdf'], 'Resolution', 300) % Use vector for best quality    


%% ** Prediction Maps 
close all; clear; clc

% CONFIGURATION
year_of_interest = 2020;
selected_regions = {'Alaska', 'Russia'}; % You can add more regions
horizons = [0, 1, 2, 5];
bin_edges = [0, 60, 70, 80, 90, 100, inf];
bin_labels = {'<60', '60-70', '70-80', '80-90', '90-100', '>100'};
cmap = parula(numel(bin_labels)); % Choose a nice colormap

% LOOP OVER HORIZONS
for h = horizons
    % Load test prediction file
    test_file = sprintf('ALT_hrz%d_TestPredictions.csv', h);
    data = readtable(test_file);

    % Filter by year and region
    mask = data.Year == year_of_interest & ismember(data.MainRegion, selected_regions);
    subset = data(mask, :);
    
    if isempty(subset)
        fprintf('⚠️ No data for horizon %d at year %d\n', h, year_of_interest);
        continue
    end

    % Bin the predicted ALT values
    bin_idx = discretize(subset.y_tepred_ensemble, bin_edges);

    % Create figure
    figure('Name', sprintf('ALT Predictions - Year %d - Horizon %d', year_of_interest, h), ...
        'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.05 0.05 0.9 0.85]);

    % Plot basemap
    ax = axesm('MapProjection','eqdcylin','MapLatLimit',[30 80],'MapLonLimit',[-180 180]);
    setm(ax, 'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on');
    geoshow('landareas.shp', 'FaceColor', [0.95 0.95 0.95]);

    % Scatter ALT values
    hold on
    for i = 1:numel(bin_labels)
        idx = bin_idx == i;
        scatterm(subset.Latitude(idx), subset.Longitude(idx), 20, ...
            'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', 'DisplayName', bin_labels{i});
    end

    % Final touches
    title(sprintf('Predicted ALT at Year %d (Horizon +%d years)', year_of_interest, h), ...
        'FontWeight', 'bold', 'FontSize', 14);
    legend('Location', 'eastoutside');
end

%% ** Focused Prediction maps

close all; clear; clc

% CONFIGURATION
year_of_interest = 2020;
selected_subregions = {'NorthSlope', 'WestSiberia'};  % Focused zoom regions
horizons = [0, 1, 2, 5];
bin_edges = [0, 60, 70, 80, 90, 100, inf];
bin_labels = {'<60', '60-70', '70-80', '80-90', '90-100', '>100'};
cmap = parula(numel(bin_labels));

% TIGHTER REGION BOUNDARIES (based on actual ALT coverage zones)
subregion_bounds = struct( ...
    'NorthSlope',    [68 71; -156 -146], ...  % Deadhorse, Toolik area
    'WestSiberia',   [60 68; 65 80]); %...        % Example zone with MODIS coverage
    % Add more: e.g., 'EasternSiberia', 'InteriorAlaska', 'Chukotka'


% LOOP OVER HORIZONS
for h = horizons
    % Load test prediction file
    test_file = sprintf('ALT_hrz%d_TestPredictions.csv', h);
    data = readtable(test_file);

    % LOOP OVER SELECTED SUBREGIONS
    for r = 1:numel(selected_subregions)
        name = selected_subregions{r};
        bounds = subregion_bounds.(name);
        lat_lim = bounds(1, :);
        lon_lim = bounds(2, :);

        % Filter by year and spatial box
        mask = data.Year == year_of_interest & ...
               data.Latitude >= lat_lim(1) & data.Latitude <= lat_lim(2) & ...
               data.Longitude >= lon_lim(1) & data.Longitude <= lon_lim(2);
        subset = data(mask, :);

        if isempty(subset)
            fprintf('⚠️ No data in %s at horizon %d, year %d\n', name, h, year_of_interest);
            continue
        end

        % Bin the predicted ALT values
        bin_idx = discretize(subset.y_tepred_ensemble, bin_edges);

        % Create figure
        figure('Name', sprintf('%s - Year %d - Horizon %d', name, year_of_interest, h), ...
            'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.8 0.85]);

        % Plot zoomed map
        ax = axesm('eqdcylin', 'MapLatLimit', lat_lim, 'MapLonLimit', lon_lim);
        setm(ax, 'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on');
        geoshow('landareas.shp', 'FaceColor', [0.95 0.95 0.95]);

        % Plot ALT predictions
        hold on
        for i = 1:numel(bin_labels)
            idx = bin_idx == i;
            scatterm(subset.Latitude(idx), subset.Longitude(idx), 25, ...
                'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', 'DisplayName', bin_labels{i});
        end

        % Final touches
        title(sprintf('Predicted ALT in %s – Year %d (Horizon +%d)', name, year_of_interest, h), ...
              'FontWeight', 'bold', 'FontSize', 14);
        legend('Location', 'eastoutside');
    end
end


%% Functions
function categoryTable = assignEnvironmentalCategories(latitudes, longitudes)
    % Load the permafrost shapefile
    shapefile = 'ggd318_map_circumarctic\permaice.shp'; 
    permafrost = shaperead(shapefile, 'UseGeoCoords', false);
    shape_info = shapeinfo(shapefile);
    source_proj = shape_info.CoordinateReferenceSystem;

    % Project the polygon coordinates
    [X, Y] = deal({permafrost.X}, {permafrost.Y});
    [polyLat, polyLon] = cellfun(@(x, y) projinv(source_proj, x, y), X, Y, 'UniformOutput', false);
    combo_values = {permafrost.COMBO};

    % Initialize results
    n = length(latitudes);
    permafrostType = repmat("Unknown", n, 1);
    groundIceType = repmat("Unknown", n, 1);
    landformType = repmat("Unknown", n, 1);

    % Loop through each point
    for i = 1:n
        lat = latitudes(i);
        lon = longitudes(i);
        found = false;

        % Check which polygon it falls into
        for j = 1:length(polyLat)
            if inpolygon(lon, lat, polyLon{j}, polyLat{j})
                code = lower(combo_values{j});
                if numel(code) >= 3
                    % Extract and label types
                    switch code(1)
                        case 'c', permafrostType(i) = "Continuous";
                        case 'd', permafrostType(i) = "Discontinuous";
                        case 's', permafrostType(i) = "Sporadic";
                        case 'i', permafrostType(i) = "Isolated";
                    end
                    switch code(2)
                        case 'h', groundIceType(i) = "High";
                        case 'm', groundIceType(i) = "Medium";
                        case 'l', groundIceType(i) = "Low";
                    end
                    switch code(3)
                        case 'f', landformType(i) = "Lowlands";
                        case 'r', landformType(i) = "Mountains";
                    end
                end
                found = true;
                break;
            end
        end

        if ~found
            % Optional: log unmatched coordinates
            fprintf("No match for site at (%.2f, %.2f)\n", lat, lon);
        end
    end

    % Create table for output
    categoryTable = table(latitudes, longitudes, permafrostType, groundIceType, landformType);
end


function corrs = chatterjee_corr_vectorized(X, y)
% CHATTERJEE_CORR_VECTORIZED Compute Chatterjee correlation between each column in X and y
%   corrs = chatterjee_corr_vectorized(X, y)
%   Inputs:
%       X - matrix of predictors (n x p)
%       y - target vector (n x 1)
%   Output:
%       corrs - 1 x p vector of Chatterjee coefficients

    % Validate input
    if istable(X)
        X = table2array(X);
    end

    [n, p] = size(X);
    y = y(:);

    if length(y) ~= n
        error('X and y must have the same number of rows.');
    end

    % Remove rows with NaNs
    rowsToUse = all(~isnan(X), 2) & ~isnan(y);
    X = X(rowsToUse, :);
    y = y(rowsToUse);
    n = size(X, 1);

    % Rank y once
    yranks = tiedrank(y);
    
    % Initialize result
    corrs = zeros(1, p);

    for j = 1:p
        xj = X(:, j);

        % Sort xj and get sorted yranks
        [~, sortIdx] = sort(xj);
        yr_sorted = yranks(sortIdx);

        % Compute absolute differences
        diffs = abs(diff(yr_sorted));
        corrs(j) = 1 - 3 * sum(diffs) / (n^2 - 1);
    end
end

% ---------- Metric Functions ----------

function val = rmse(y, yhat)
    val = sqrt(mean((y - yhat).^2));
end

function val = mae(y, yhat)
    val = mean(abs(y - yhat));
end

function val = r2score(y, yhat)
    ss_res = sum((y - yhat).^2);
    ss_tot = sum((y - mean(y)).^2);
    val = 1 - ss_res / ss_tot;
end

% ---------- Confidence Interval ----------

function [mean_metric, ci] = compute_ci(y, yhat, metric_fn)
    rng(42);  % Set once before bootstrapping
    n_boot = 1000;
    n = length(y);
    vals = zeros(n_boot, 1);
    for i = 1:n_boot
        idx = randi(n, n, 1);
        vals(i) = metric_fn(y(idx), yhat(idx));
    end
    mean_metric = mean(vals);
    ci = prctile(vals, [2.5, 97.5]);
end

% ---------- Bar Plotting Function ----------

function plot_performance_by_horizon(tbl, methods)
% Colors
colors = [
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.8500 0.3250 0.0980;
    0.4940 0.1840 0.5560;
    0.9290 0.6940 0.1250;
    0.6350 0.0780 0.1840;
];

tbl.Method = categorical(tbl.Method, methods, 'Ordinal', true);
horizons = unique(tbl.Horizon);
hz_labels = compose('+%d days', horizons);
metrics = {'Test_RMSE', 'Test_MAE', 'Test_R2'};
metric_labels = {'RMSE', 'MAE', 'R$^2$'};

figure('Units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(3, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for k = 1:length(metrics)
    nexttile;
    hold on;

    data_vals = NaN(length(horizons), length(methods));
    data_err_low = NaN(length(horizons), length(methods));
    data_err_up = NaN(length(horizons), length(methods));

    for h = 1:length(horizons)
        for m = 1:length(methods)
            row = strcmp(string(tbl.Method), methods{m}) & tbl.Horizon == horizons(h);
            val = tbl{row, metrics{k}};
            lwr = tbl{row, [metrics{k}, '_L']};
            upr = tbl{row, [metrics{k}, '_U']};

            data_vals(h, m) = val;
            data_err_low(h, m) = lwr;
            data_err_up(h, m) = upr;
        end
    end

    b = bar(data_vals, 'grouped');
    for m = 1:length(methods)
        b(m).FaceColor = colors(m, :);
        b(m).EdgeColor = 'none';
    end

    % Add asymmetric error bars
    [ngroups, nbars] = size(data_vals);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, data_vals(:,i), (data_err_up(:,i) - data_err_low(:,i))/2, ...
            'k', 'linestyle', 'none', 'LineWidth', 1.2, 'CapSize', 10);

        % Add latex formatted values
        for j = 1:ngroups
            val = data_vals(j, i);
            lwr = data_err_low(j, i);
            upr = data_err_up(j, i);
            % label = sprintf('$%.2f^{%.2f}_{%.2f}$', val, upr, lwr);
            label = sprintf('%.2f±%.2f', val, (upr-lwr)/2);

            if i ~= 3 || (k == 3 && i==3) %val > 1.5 * (upr - lwr)
                yloc = val / 2; % inside
                color = 'w';
            else
                if k ~= 1
                    yloc = val + data_err_up(j, i) + .02 * max(data_vals(:));
                    color = 'k';
                    x(j) = x(j);
                else
                    yloc = val / 2;%val + data_err_up(j, i) + 0.01 * max(data_vals(:));
                    color = 'w';
                    x(j) = x(j) +.03;
                end
            end

            text(x(j), yloc, label, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 15, 'Interpreter', 'latex', ...
                'Color', color, Rotation=90);
        end
    end

    set(gca, 'XTick', 1:length(hz_labels), 'XTickLabel', [], ...
        'FontSize', 16, 'TickLabelInterpreter','latex');
    ylabel(metric_labels{k}, 'Interpreter', 'latex', 'FontSize', 20);
    grid on; box on;

    if k == 3
        ylim([0 1]);
        set(gca, 'XTick', 1:length(hz_labels), 'XTickLabel', hz_labels, ...
            'FontSize', 16, 'TickLabelInterpreter','latex');
    end
end

lgd = legend(methods, 'Location', 'northoutside', 'Orientation', 'horizontal');
lgd.Layout.Tile = 'north';
set(lgd, 'Interpreter', 'latex', 'FontSize', 16);
end


function [corr_cat, coordinates, year, LC, topography, vegetation, PFG, SFG, NFG, ERA5] = computeCategoryMeans(xi_corr,bias)
    % Function to compute category means for a given xi_corr variable.
    % xi_corr: Input correlation matrix.
    % bias: could be 3 if the xi_corr does not contain the first three
    % categorical features (mainRegion, etc..)
    coordinates = xi_corr(:, 1+bias:2+bias);
    year = xi_corr(:, 3+bias);
    LC = xi_corr(:, 4+bias:8+bias);
    topography = xi_corr(:, 9+bias:11+bias);
    vegetation = xi_corr(:, 12+bias:13+bias);
    PFG = xi_corr(:, [14+bias, 17+bias:32+bias]);
    SFG = xi_corr(:, [15+bias, 33+bias:48+bias]);
    NFG = xi_corr(:, [16+bias, 49+bias:64+bias]);
    ERA5 = xi_corr(:, 65+bias:79+bias);

    corr_cat1 = mean(xi_corr(:, 1+bias:2+bias), 1); % coordinates
    corr_cat2 = mean(xi_corr(:, 3+bias), 1);   % year
    corr_cat3 = mean(xi_corr(:, 4+bias:8+bias), 1); % LC
    corr_cat4 = mean(xi_corr(:, 9+bias:11+bias), 1); % topography
    corr_cat5 = mean(xi_corr(:, 12+bias:13+bias), 1); % vegetation
    corr_cat6 = mean(xi_corr(:, [14+bias, 17+bias:32+bias]), 1); % PFG
    corr_cat7 = mean(xi_corr(:, [15+bias, 33+bias:48+bias]), 1); % SFG
    corr_cat8 = mean(xi_corr(:, [16+bias, 49+bias:64+bias]), 1); % NFG
    corr_cat9 = mean(xi_corr(:, 65+bias:79+bias), 1); % ERA5

    % Combine results into final output
    corr_cat = [corr_cat1 corr_cat2 corr_cat3 corr_cat4 ...
                corr_cat5 corr_cat6 corr_cat7 corr_cat8, corr_cat9];
end


% Helper function: Compute Entropy
function ent = compute_entropy(seq)
    seq = seq(~isnan(seq)); % Remove NaNs
    if isempty(seq)
        ent = NaN;
        return;
    end
    prob_dist = histcounts(seq, 'Normalization', 'probability');
    prob_dist = prob_dist(prob_dist > 0); % Remove zero probabilities
    ent = -sum(prob_dist .* log2(prob_dist));
end


% Helper function: Compute Lempel-Ziv Complexity
function lzc = lempel_ziv_complexity(seq)
    seq = normalize(seq, 'range', [0, 255]);
    bin_seq = dec2bin(uint8(seq));
    unique_subs = {}; idx = 1;
    while idx <= length(bin_seq)
        for j = idx:length(bin_seq)
            sub_seq = bin_seq(idx:j);
            if ~any(strcmp(unique_subs, sub_seq))
                unique_subs{end+1} = sub_seq;
                break;
            end
        end
        idx = j + 1;
    end
    lzc = length(unique_subs) / length(bin_seq);
end


