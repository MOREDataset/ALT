%% **Introduction Permafrost types CALM MAPs
clear
clc
% Load the data file
filePath = 'CALM_ALT_all.xlsx';

% Read the data starting from the correct header row
opts = detectImportOptions(filePath, 'Sheet', 'Sheet1');
data = readtable(filePath, opts);

% Extract relevant columns
locations = data.Location;
years = data.Year;
maxALT = data.Max;

% Remove rows with missing data
validIdx = ~iALT_DatasetDistributions_V3smissing(locations) & ~ismissing(years) & ~ismissing(maxALT);
data = data(validIdx, :);
CALMlocations = unique([data.Lat data.Long],"rows");

%*Polar projection permafrost distribution map with CALM sites: Load Permafrost Shapefile and Convert Projection
shapefile = 'ggd318_map_circumarctic\permaice.shp'; 
permafrost = shaperead(shapefile, 'UseGeoCoords', false); % Faster reading

% Extract projection information
shape_info = shapeinfo(shapefile);
source_proj = shape_info.CoordinateReferenceSystem;

% Convert all projected coordinates to lat/lon **in one step**
[X, Y] = deal({permafrost.X}, {permafrost.Y});
[lat, lon] = cellfun(@(x, y) projinv(source_proj, x, y), X, Y, 'UniformOutput', false);
CALMlocations = unique([data.Lat data.Long],"rows");

%*Categorize Data: Permafrost, Ground Ice, Landform**
categories = struct( ...
    'permafrost', struct( ...
        'continuous', {'c'}, ...
        'discontinuous', {'d'}, ...
        'sporadic', {'s'}, ...
        'isolated', {'i'}), ...
    'ground_ice', struct( ...
        'high', {'h'}, ...
        'medium', {'m'}, ...
        'low', {'l'}), ...
    'landform', struct( ...
        'lowlands', {'f'}, ...
        'mountains', {'r'}));

% Define Colors for Each Layer
colors.permafrost = struct( ...
    'continuous', [0.0,0.2,0.8], ...
    'discontinuous', [0.2,0.4,0.9], ...
    'sporadic', [0.4,0.6,1.0], ...
    'isolated', [0.7,0.85,1.0]);

colors.ground_ice = struct( ...
    'high', [0.0,0.4784,0.4667], ...       % Dark Green
    'medium', [0.2,0.7,0.7], ...     % Green
    'low', [0.6,0.85,0.85]);           % Light Green

colors.landform = struct( ...
    'lowlands', [0.0,0.5,0.0], ...  % Brown
    'mountains', [198,130,68]./255);  % Gray

% Extract all `COMBO` values at once
combo_values = {permafrost.COMBO};

% Initialize category assignments
num_polygons = length(permafrost);
permafrost_types = repmat({'other'}, num_polygons, 1);
ground_ice_types = repmat({'other'}, num_polygons, 1);
landform_types = repmat({'other'}, num_polygons, 1);

% Assign categories
for i = 1:num_polygons
    combo_code = lower(combo_values{i}); % Convert to lowercase

    % Extract permafrost type
    if numel(combo_code) > 2 
        permafrost_types{i} = combo_code(1);
        ground_ice_types{i} = combo_code(2);
        landform_types{i} = combo_code(3);
    end
end

%*Plot the Permafrost Extent Map**
figure('units','normalized','OuterPosition',[0 0 1 1]);
ax = axesm('stereo', 'MapLatLimit', [50 90], 'Frame', 'on', 'Grid', 'on');
setm(ax, 'MeridianLabel', 'off', 'ParallelLabel', 'on');

% Load land shapefile (Natural Earth land data)
land = shaperead('landareas.shp', 'UseGeoCoords', true);

% Plot land areas in light gray
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'k');

field = fieldnames(categories.permafrost)';
for i = 1:numel(field)
    category =  categories.permafrost.(field{i});
    idx = strcmp(permafrost_types, category);
    if any(idx)
        h1(i) = geoshow([lat{idx}], [lon{idx}], 'DisplayType', 'polygon', ...
            'FaceColor', colors.permafrost.(field{i}), 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    end
end
lgd1 = legend([h1(1) h1(2) h1(3) h1(4)],{'Continuous (90-100\%)', 'Discontinuous (50-90\%)', 'Sporadic (10-50\%)', 'Isolated (0-10\%)'});
lgd1.Title.String = 'Permafrost Distribution';
set(lgd1,'Location', 'northeastoutside', 'Orientation', 'vertical', 'Interpreter','latex','FontSize',18);

% Plot CALM sites
geoshow(CALMlocations(:,1), CALMlocations(:,2), 'DisplayType', 'point', 'Marker', 'pentagram', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0.6350 0.0780 0.1840],'MarkerSize',10, 'DisplayName','CALM sites');

% exportgraphics(gcf, ['C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\','CALM_permaMap1.pdf'],'Resolution',300)

%*Plot Ground Ice Content Map**
figure('units','normalized','OuterPosition',[0 0 1 1]);
ax = axesm('stereo', 'MapLatLimit', [50 90], 'Frame', 'on', 'Grid', 'on');
setm(ax, 'MeridianLabel', 'off', 'ParallelLabel', 'on');

% Load land shapefile (Natural Earth land data)
land = shaperead('landareas.shp', 'UseGeoCoords', true);

% Plot land areas in light gray
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'k');

field = fieldnames(categories.ground_ice)';
for i = 1:numel(field)
    category =  categories.ground_ice.(field{i});
    idx = strcmp(ground_ice_types, category);
    if any(idx)
        h2(i) = geoshow([lat{idx}], [lon{idx}], 'DisplayType', 'polygon', ...
            'FaceColor', colors.ground_ice.(field{i}), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
end
lgd2 = legend([h2(1) h2(2) h2(3)],{'High ($>$20\%)', 'Medium (10-20\%)', 'Low (0-10\%)'});
lgd2.Title.String = 'Ground Ice Content';
set(lgd2,'Interpreter','latex','FontSize',14, 'Location','northeastoutside','Orientation','vertical');

% Plot CALM sites
geoshow(CALMlocations(:,1), CALMlocations(:,2), 'DisplayType', 'point', 'Marker', 'pentagram', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0.6350 0.0780 0.1840],'MarkerSize',10);%, 'DisplayName','CALM sites'

% exportgraphics(gcf, ['C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\','CALM_permaMap2_noMRKcaption.pdf'],'Resolution',300)

%**Plot Landform**
figure('units','normalized','OuterPosition',[0 0 1 1]);
ax = axesm('stereo', 'MapLatLimit', [50 90], 'Frame', 'on', 'Grid', 'on');
setm(ax, 'MeridianLabel', 'off', 'ParallelLabel', 'on');

% Load land shapefile (Natural Earth land data)
land = shaperead('landareas.shp', 'UseGeoCoords', true);

% Plot land areas in light gray
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'k');

field = fieldnames(categories.landform)';
for i = 1:numel(field)
    category =  categories.landform.(field{i});
    idx = strcmp(landform_types, category);
    if any(idx)
        h3(i) = geoshow([lat{idx}], [lon{idx}], 'DisplayType', 'polygon', ...
            'FaceColor', colors.landform.(field{i}), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end
end

% Add title & legend
% title('Circum-Arctic Permafrost Distribution: Extent, Ground Ice, and Landform');

lgd3 = legend([h3(1) h3(2)], {sprintf('Thick overburden in lowlands\n\\& depressions'),sprintf('Thin overburden in mountains\n\\& plateaus')},'Location', 'northeastoutside', 'Orientation', 'vertical', 'Interpreter','latex','FontSize',18);
lgd3.Title.String = 'Landform types';
set(lgd3,'Interpreter','latex','FontSize',14);

% Plot CALM sites
geoshow(CALMlocations(:,1), CALMlocations(:,2), 'DisplayType', 'point', 'Marker', 'pentagram', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0.6350 0.0780 0.1840],'MarkerSize',10);%, 'DisplayName','CALM sites'

% exportgraphics(gcf, ['C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\','CALM_permaMap3_noMRKcaption.pdf'],'Resolution',300)

disp('Permafrost, Ground Ice, and Landform mapped successfully.');

%% 1** CALM data with Standardized Location Names, Regularized ALT data for regression
clc; clear;
filename = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\CALM_ALT_all.xlsx';
data = readtable(filename);
fprintf('Loaded the CALM dataset\n');

% Step 1: Standardize Location Names (Merge sites with same coordinates)
[~, ia, ic] = unique([data.Lat, data.Long], "rows"); % Find unique coordinates
standardized_names = data.Location(ia); % Pick one name per coordinate
data.Location = standardized_names(ic); % Assign standardized names
fprintf('Standardized location names for sites with duplicate coordinates.\n');

% Add regions to the dataset
regions = {'Svalbard','North East Siberia','North East Siberia','Canada','North East Siberia',...
           'North East Siberia','Alaska Interior','North East Siberia','Alaska Interior','Alaska North Slope',...
           'Russian European north','Alaska North Slope','Canada','Russian European north','Alaska Interior',...
           'Canada','Central Siberia','Central Siberia','Canada','North East Siberia','North East Siberia',...
           'Chukotka','Alaska North Slope','Central Siberia','Central Siberia','Alaska Interior',...
           'North East Siberia','Alaska Interior','Alaska Seward Peninsula','Canada','Alaska North Slope',...
           'Alaska North Slope','Alaska Interior','Svalbard','Alaska Interior','West Siberia','Alaska Interior',...
           'Alaska Interior','Alaska North Slope','Alaska North Slope','Alaska Interior','West Siberia',...
           'West Siberia','West Siberia','West Siberia','West Siberia','Alaska North Slope','West Siberia',...
           'Alaska Interior','Central Siberia','Alaska North Slope','Alaska North Slope','North East Siberia',...
           'Russian European north','North East Siberia','Alaska Interior','Alaska Seward Peninsula','Chukotka',...
           'West Siberia','North East Siberia','North East Siberia','Chukotka','Chukotka','Alaska North Slope',...
           'North East Siberia','West Siberia','Alaska Interior','Alaska Interior','Chukotka','North East Siberia',...
           'North East Siberia','North East Siberia','North East Siberia','North East Siberia','West Siberia',...
           'Central Siberia','Alaska Interior','Alaska Interior','West Siberia','West Siberia','Alaska Interior',...
           'North East Siberia','Canada','West Siberia','Alaska Interior','Kamchatka','Kamchatka','Kamchatka',...
           'Alaska Interior','Alaska North Slope','West Siberia','North East Siberia','Russian European north',...
           'Alaska Interior','Central Siberia','Russian European north','Alaska Interior','Canada','Canada',...
           'Alaska North Slope','Central Siberia','Svalbard','Canada','Russian European north','West Siberia',...
           'West Siberia','West Siberia','West Siberia','West Siberia','West Siberia','West Siberia',...
           'Alaska North Slope','Alaska Interior','Svalbard','Svalbard'};
      

locsRegions = table(unique(data.Location), regions','VariableNames',{'Location','Region'});
data = join(data,locsRegions,"Keys","Location");
data = movevars(data, "Region", "Before", "Location");

%%corrections
% tmp = data(strcmp(data.Region,'Svalbard'),:);
% [~,ia,~] = unique([tmp.Lat tmp.Long],'rows')
% tmp(ia,:)
data.Region(contains(data.Location,'Sweden')) = cellstr(repmat('Sweden', length(data.Region(contains(data.Location,'Sweden'))), 1));
data.Region(contains(data.Location,'Zackenberg')) = cellstr(repmat('Greenland', length(data.Region(contains(data.Location,'Zackenberg'))), 1));

% Define main categories and their respective sub-regions
main_categories = {'Alaska', 'Canada', 'Russia', 'Svalbard', 'Greenland', 'Sweden'};
category_mapping = containers.Map(...
    {'Alaska North Slope', 'Alaska Interior', 'Alaska Seward Peninsula',...
     'Central Siberia', 'West Siberia', 'North East Siberia', 'Russian European north', 'Chukotka', 'Kamchatka',...
     'Canada', 'Svalbard', 'Greenland', 'Sweden'}, ...
    {'Alaska', 'Alaska', 'Alaska',...
     'Russia', 'Russia', 'Russia', 'Russia', 'Russia', 'Russia',...
     'Canada', 'Svalbard', 'Greenland', 'Sweden'});

categories_colors = lines(length(main_categories)); % Assign unique colors

% Assign main categories based on sub-region
data.MainRegion = cell(height(data), 1);
for i = 1:height(data)
    if isKey(category_mapping, data.Region{i})
        data.MainRegion{i} = category_mapping(data.Region{i});
    else
        data.MainRegion{i} = 'Other';
    end
end

% Column names
site_column = 'Location';
year_column = 'Year';
alt_column = 'Max';

% Extract unique locations based on lat/lon
[~, ia, ~] = unique([data.Lat, data.Long], "rows");
sites = data.Location(ia);
unique_lats = data.Lat(ia);
unique_lons = data.Long(ia);

% Initialize a table to store all preprocessed sequences
preprocessed_data = table();

% Initialize storage for results
result_sites = {};
start_years = [];
end_years = [];
lon = [];
lat = [];
region = [];
mainRegion = [];

% Iterate over each site (unique lat/lon)
for i = 1:length(sites)
    fprintf('*%s:\n',sites{i})
    % Extract rows for this specific location
    site_rows = (data.Lat == unique_lats(i)) & (data.Long == unique_lons(i));
    site_data = data(site_rows, :);

    % Sort data by year
    site_years = site_data.Year;
    site_alt = site_data.Max;
    [sorted_years, sort_idx] = sort(site_years);
    sorted_alt = site_alt(sort_idx);

    % Define full range of years for this specific location
    min_year = min(sorted_years);
    max_year = max(sorted_years);
    full_years = (min_year:max_year)'; % Unique full range per site

    % Regularize ALT time series (align to full year range)
    regularized_alt = NaN(length(full_years), 1);
    [~, loc] = ismember(sorted_years, full_years);
    regularized_alt(loc) = sorted_alt; % Assign existing values

    % Interpolate missing years using pchip
    regularized_alt = fillmissing(regularized_alt, 'pchip');

    % Create a table row for this site
    site_table = table(repmat(site_data.MainRegion(1), length(full_years), 1),...
        repmat(site_data.Region(1), length(full_years), 1),...
        repmat(site_data.Location(1), length(full_years), 1),...
        repmat(site_data.Lat(1), length(full_years), 1), ...
        repmat(site_data.Long(1), length(full_years), 1), ...
        full_years, regularized_alt, ...
        'VariableNames', {'MainRegion','Region','Location','Latitude', 'Longitude', 'Year', 'ALT_Max'});

    % Append to the main table
    preprocessed_data = [preprocessed_data; site_table];
end

writetable(preprocessed_data,'ALT_preprocessedDataset.csv');

fprintf('Standardized locations and regularized ALT time series (per location).\n');

%% ** ALT Dataset Analysis: Range of Years, Samples, Unique Locations, and Statistical Metrics
clear
clc

% Load the dataset
filename = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_preprocessedDataset.csv';
data = readtable(filename);
data.MainRegion = categorical(data.MainRegion);
mainregions = unique(data.MainRegion);

% Preallocate results table
final_table = table( ...
    categorical(), ... % MainRegion
    zeros(0,1), ...    % num_samples
    zeros(0,1), ...    % unique_lat_lon
    strings(0,1), ...  % minmaxYears
    zeros(0,1), ...    % maxrange_years
    zeros(0,1), ...    % minrange_years
    zeros(0,1), ...    % mean_ALT
    zeros(0,1), ...    % std_ALT
    zeros(0,1), ...    % trend_strength
    zeros(0,1), ...    % seasonality_strength
    zeros(0,1), ...    % lzc_ALT
    'VariableNames', {'MainRegion','num_samples','unique_lat_lon','minmaxYears', ...
                      'maxrange_years','minrange_years','mean_ALT','std_ALT', ...
                      'trend_strength','seasonality_strength','lzc_ALT'});

% Loop over main regions
for n = 1:length(mainregions)
    region_name = mainregions(n);
    temp = data(data.MainRegion == region_name, :);
    lat_lon_groups = findgroups(temp.Latitude, temp.Longitude);

    % Compute per-site statistics
    results = table();
    results.range_years = splitapply(@(y) length(unique(y)), temp.Year, lat_lon_groups);
    results.mean_ALT = splitapply(@mean, temp.ALT_Max, lat_lon_groups);
    results.std_ALT = splitapply(@std, temp.ALT_Max, lat_lon_groups);
    results.trend_strength = splitapply(@mann_kendall_trend, temp.ALT_Max, lat_lon_groups);
    results.seasonality_strength = splitapply(@seasonality_strength_safe, temp.ALT_Max, lat_lon_groups);
    results.lzc_ALT = splitapply(@lempel_ziv_complexity, temp.ALT_Max, lat_lon_groups);

    % Aggregate statistics
    final_table.MainRegion(n)            = region_name;
    final_table.num_samples(n)           = height(temp);
    final_table.unique_lat_lon(n)        = numel(unique(lat_lon_groups));
    final_table.minmaxYears(n)           = sprintf('%d-%d', min(temp.Year), max(temp.Year));
    final_table.maxrange_years(n)        = max(results.range_years);
    final_table.minrange_years(n)        = min(results.range_years);
    final_table.mean_ALT(n)              = mean(results.mean_ALT, 'omitnan');
    final_table.std_ALT(n)               = mean(results.std_ALT, 'omitnan');
    final_table.trend_strength(n)        = mean(results.trend_strength, 'omitnan');
    final_table.seasonality_strength(n)  = mean(results.seasonality_strength, 'omitnan');
    final_table.lzc_ALT(n)               = mean(results.lzc_ALT, 'omitnan');
end

% Display results
disp('Final Averaged ALT Dataset Statistics by Main Region:');
disp(final_table);


%% ** Train and Test Splitting: Temporal Test Set (80/20 Split)
clear;
clc;

for hrz = [0, 1, 2, 5]
    file = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_preprocessedDataset.csv';
    opts = detectImportOptions(file);
    data = readtable(file, opts);
    data.Year = [];

    file = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_preprocessedDataset_wTopography.csv';
    opts = detectImportOptions(file);
    dataTopography = readtable(file, opts);
    dataTopography(:,["Elevation", "Slope", "Aspect"]) = fillmissing(dataTopography(:,["Elevation", "Slope", "Aspect"]),"constant",0);
    fprintf(', Topography data size %d',height(dataTopography))

    file = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_preprocessedDataset_wLC.csv';
    opts = detectImportOptions(file);
    dataLC = readtable(file, opts);
    fprintf(', LC data size %d',height(dataLC))

    file = ['C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_preprocessedDataset_wVegetation_', num2str(hrz),'Years.csv'];
    opts = detectImportOptions(file);
    dataVegetation_Xyears = readtable(file, opts);
    fprintf(', Vegetation data size %d',height(dataVegetation_Xyears))

    % file = ['C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_preprocessedDataset_wSnow_', num2str(hrz),'Years.csv'];
    % opts = detectImportOptions(file);
    % dataSnow_Xyears = readtable(file, opts);
    % fprintf(', Snow data size %d',height(dataSnow_Xyears))

    file = ['C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_preprocessedDataset_wERA5_', num2str(hrz),'Years.csv'];
    opts = detectImportOptions(file);
    dataERA5_Xyears = readtable(file, opts);
    dataERA5_Xyears.Latitude = [];
    dataERA5_Xyears.Longitude= [];
    dataERA5_Xyears(:,6:end) = fillmissing(dataERA5_Xyears(:,6:end), 'constant', 0);
    fprintf(', ERA5 data size %d',height(dataERA5_Xyears))

    file = ['C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\CALM\ALT_preprocessedDataset_wFG_', num2str(hrz),'Years.csv'];
    opts = detectImportOptions(file);
    dataFG_Xyears = readtable(file, opts);
    fprintf(', FG data size %d',height(dataFG_Xyears))
    dataFG_Xyears.Lat = [];
    dataFG_Xyears.Lon= [];
    dataFG_Xyears(:,11:end) = fillmissing(dataFG_Xyears(:,11:end), 'constant', 0);

    data = [data, dataVegetation_Xyears(:,6), dataLC(:,8:end), dataTopography(:,8:end), dataVegetation_Xyears(:,8:end), dataFG_Xyears(:,8:end), dataERA5_Xyears(:,8:end)];% dataERA5_Xyears(:,8:end),

    % Rename variables if necessary
    data.Properties.VariableNames{4} = 'Lat';
    data.Properties.VariableNames{5} = 'Long';

    % Step 1: Find unique locations
    uniqueLocations = unique(data(:, {'Lat', 'Long'}), 'rows');

    % Compute number of records per location
    locationCounts = varfun(@height, data, 'GroupingVariables', {'Lat', 'Long'});
    locationCounts = sortrows(locationCounts, 'GroupCount', 'descend');

    % Step 2: Determine 80/20 split
    totalRecords = height(data);
    numTestRecords = round(0.2 * totalRecords);  % 20% for testing
    remainingData = sortrows(data, 'Year', 'descend'); % Sort by year (most recent first)

    % Step 3: Select 20% of data for testing while considering locations
    temporalTestData = [];
    trainData = remainingData;
    uniqueLocs = unique(remainingData(:, {'Lat', 'Long'}), 'rows');

    for i = 1:height(uniqueLocs)
        lat = uniqueLocs.Lat(i);
        lon = uniqueLocs.Long(i);

        % Extract all records for this location
        locData = remainingData(remainingData.Lat == lat & remainingData.Long == lon, :);

        % Select recent records per location while keeping a total 20% split
        numLocRecords = height(locData);
        numLocTestRecords = round(0.2 * numLocRecords);  % Select 20% per location

        if numLocTestRecords > 0
            testSubset = locData(1:numLocTestRecords, :); % Select the most recent records
            trainSubset = setdiff(locData, testSubset);

            % Append to temporal test and train sets
            temporalTestData = [temporalTestData; testSubset];
            trainData = setdiff(trainData, testSubset);
        end
    end

    % ✅ Ensure No Unseen Locations in Test Set
    trainLocations = unique([trainData.Lat, trainData.Long], 'rows');
    temporalTestLocations = unique([temporalTestData.Lat, temporalTestData.Long], 'rows');

    TemporalUnseen = setdiff(temporalTestLocations, trainLocations, 'rows'); % Only locations not in training set
    if ~isempty(TemporalUnseen)
        % Extract records with unseen locations
        unseenTemporalData = innerjoin(temporalTestData, array2table(TemporalUnseen, 'VariableNames', {'Lat', 'Long'}), 'Keys', {'Lat', 'Long'});

        % Remove unseen locations from the temporal test set
        temporalTestData = setdiff(temporalTestData, unseenTemporalData);

        % Add them back to the training set
        trainData = [trainData; unseenTemporalData];

        % Recompute train locations after adjustment
        trainLocations = unique([trainData.Lat, trainData.Long], 'rows');
        temporalTestLocations = unique([temporalTestData.Lat, temporalTestData.Long], 'rows');
    end

    % ✅ Compute and Display Updated Dataset Summary

    summaryTable = table(["Trainset"; "Temporal testset"], ...
        [height(trainData); height(temporalTestData)], ...
        [height(trainData)/height(data)*100; height(temporalTestData)/height(data)*100], ...
        [size(trainLocations,1); size(temporalTestLocations,1)], ...
        [size(trainLocations,1); size(intersect(trainLocations, temporalTestLocations, 'rows'),1)], ...
        [100; size(intersect(trainLocations, temporalTestLocations, 'rows'),1) / size(trainLocations,1) * 100], ...
        [0; size(setdiff(temporalTestLocations, trainLocations, 'rows'),1)], ...
        [string(min(trainData.Year)) + " - " + string(max(trainData.Year)); ...
        string(min(temporalTestData.Year)) + " - " + string(max(temporalTestData.Year))], ...
        'VariableNames', {'Dataset', 'Total_Records', 'Total_Records rate', 'Total unique locations', 'Similar_location2Train', 'Similar_location2Train rate', 'UnSimilar_location2Train','Year_Range'});

    disp(summaryTable);

    % Save datasets % make sure the ERA5 data is used for this,
    % otherwise, remove the addition in the filename
    save(['ALT_TrainData_targetswERA5_' num2str(hrz)],"trainData");
    save(['ALT_TemporalTestData_targetswERA5_' num2str(hrz)],"temporalTestData");
end



%% ** Preprocess and export train and test sets per horizon
clear
clc


for hrz = [0 1 2 5] %
    fprintf('horizon: +%d years\n',hrz)
    load(['ALT_TrainData_targetswERA5_',num2str(hrz),'.mat']);
    trainData.Properties.VariableNames{4} = 'Latitude';
    trainData.Properties.VariableNames{5} = 'Longitude';
    load(['ALT_TemporalTestData_targetswERA5_',num2str(hrz),'.mat']);
    temporalTestData.Properties.VariableNames{4} = 'Latitude';
    temporalTestData.Properties.VariableNames{5} = 'Longitude';
    fprintf('Train size %d, temptest size %d',height(trainData),height(temporalTestData))
    
    %fill missing data in EVI/NDVI and LCx
    trainData(:,"NDVI") = fillmissing(trainData(:,"NDVI"),'constant',0);
    trainData(:,"EVI") = fillmissing(trainData(:,"EVI"),'constant',0);
    temporalTestData(:,"NDVI") = fillmissing(temporalTestData(:,"NDVI"),'constant',0);
    temporalTestData(:,"EVI") = fillmissing(temporalTestData(:,"EVI"),'constant',0);
    trainData(:,"LC_Type1") = fillmissing(trainData(:,"LC_Type1"),'constant',0);
    trainData(:,"LC_Type2") = fillmissing(trainData(:,"LC_Type2"),'constant',0);
    trainData(:,"LC_Type3") = fillmissing(trainData(:,"LC_Type3"),'constant',0);
    trainData(:,"LC_Type4") = fillmissing(trainData(:,"LC_Type4"),'constant',0);
    trainData(:,"LC_Type5") = fillmissing(trainData(:,"LC_Type5"),'constant',0);
    temporalTestData(:,"LC_Type1") = fillmissing(temporalTestData(:,"LC_Type1"),'constant',0);
    temporalTestData(:,"LC_Type2") = fillmissing(temporalTestData(:,"LC_Type2"),'constant',0);
    temporalTestData(:,"LC_Type3") = fillmissing(temporalTestData(:,"LC_Type3"),'constant',0);
    temporalTestData(:,"LC_Type4") = fillmissing(temporalTestData(:,"LC_Type4"),'constant',0);
    temporalTestData(:,"LC_Type5") = fillmissing(temporalTestData(:,"LC_Type5"),'constant',0);
    Trainset_full = movevars(trainData, "ALT_Max", "After", trainData.Properties.VariableNames{end});%dataVegetation_Xyears, dataFG_Xyears, dataSnow_Xyears, dataLST_Xyears
    Testset_Temporal_full = movevars(temporalTestData, "ALT_Max", "After", trainData.Properties.VariableNames{end});
    fprintf('\n-Before removing missing: train size %d, temptest size %d',height(Trainset_full),height(Testset_Temporal_full))

    %deal with missing data?
    Trainset_full = rmmissing(Trainset_full);
    Testset_Temporal_full = rmmissing(Testset_Temporal_full);
    fprintf('\n-After removing missing: train size %d, temptest size %d\n',height(Trainset_full),height(Testset_Temporal_full))

    % Save testing and training datasets
    writetable(Testset_Temporal_full, ['TestDatawERA5_',num2str(hrz),'Years_V3.csv']);
    writetable(Trainset_full, ['TrainDatawERA5_',num2str(hrz),'Years_V3.csv']);

    clearvars -except Testset_Temporal_full Trainset_full    %deletes all variables except X in workspace

end


%% ** Compute performance of Naive Baselines (1-2)
clear
clc
load('ALT_TrainData_targets.mat');
load('ALT_TemporalTestData_targets.mat');
testData = temporalTestData;

%%Define Prediction Horizons
horizons = [0, 1, 2, 5];

%%Define Unique Locations
locations = unique(trainData(:, {'Lat', 'Long'}), 'rows');

results = [];

for h = horizons
    %%Remove the latest `h` years from the training set
    maxYearPerLoc = groupsummary(trainData, {'Lat', 'Long'}, 'max', 'Year'); % Get max year per location
    trainDataH = trainData;

    % Loop through each unique location and remove the latest `h` years
    for i = 1:height(maxYearPerLoc)
        lat_i = maxYearPerLoc.Lat(i);
        long_i = maxYearPerLoc.Long(i);
        maxYear = maxYearPerLoc.max_Year(i);

        % Remove data for this location if its year falls within the last `h` years
        trainDataH(trainDataH.Lat == lat_i & trainDataH.Long == long_i & ...
                   trainDataH.Year > (maxYear - h), :) = [];
    end

    %%Location Agnostic, Average by Year Baseline
    yearlyMeanALT = groupsummary(trainDataH, 'Year', 'mean', 'ALT_Max');
    testData2 = outerjoin(testData, yearlyMeanALT, "Type", "left", 'Keys', 'Year', 'MergeKeys', true);
    predYearAvg = testData2.mean_ALT_Max;

    % %%Time Agnostic, Average by Unique Location Baseline
    % locationMeanALT = groupsummary(trainDataH, {'Lat', 'Long'}, 'mean', 'ALT_Max');
    % testData2 = outerjoin(testData, locationMeanALT, "Type", "left", 'Keys', {'Lat', 'Long'}, 'MergeKeys', true);
    % predLocAvg = testData2.mean_ALT_Max;

    %%Combined Baseline: Average by Year for Each Unique Location
    yearlyLocMeanALT = groupsummary(trainDataH, {'MainRegion'}, 'mean', 'ALT_Max');
    testData2 = outerjoin(testData, yearlyLocMeanALT, "Type", "left", 'Keys', {'MainRegion'}, 'MergeKeys', true);
    predYearLocAvg = testData2.mean_ALT_Max;

    %%Compute Metrics
    actual = testData.ALT_Max;
    models = {'YearAvg','RegionAvg'};%'LocAvg', };
    preds = {predYearAvg, predYearLocAvg, predLocAvg};
    
    for i = 1:length(models)
        pred = preds{i};
        validIdx = ~isnan(pred) & ~isnan(actual);
        pred = pred(validIdx);
        actualValid = actual(validIdx);

        % Metrics
        RMSE = sqrt(mean((actualValid - pred).^2));
        MAE = mean(abs(actualValid - pred));
        % MAPE = mean(abs((actualValid - pred) ./ actualValid)) * 100;
        R2 = 1 - sum((actualValid - pred).^2) / sum((actualValid - mean(actualValid)).^2);
        
        % Store Results
        results = [results; {h, models{i}, RMSE, MAE, MAPE, R2}];
    end
end

%%Convert to Table and Display
resultsTable = cell2table(results, 'VariableNames', {'Horizon', 'Model', 'RMSE', 'MAE', 'R2'});
disp(resultsTable);

%% ** train/test sets summary table
clear; clc;

horizons = [0, 1, 2, 5];  % Define prediction horizons
targetVar = 'ALT_Max';    % Change to match your target column
datetimeFormat = 'yyyy-MM-dd';  % Adjust if your dates have a different format

summaryAll = [];

for h = horizons
    fprintf('\n===== Horizon: +%d Years =====\n', h);

    trainFile = strjoin(["TrainDatawERA5_", num2str(h), "Years_V3.csv"],'');  

    % Read tables
    trainData = readtable(trainFile);

    % Summary for Train
    summaryTrain = getSummary(trainData, targetVar);

    % Combine and store
    T = table;
    T.Dataset = {'Train'};
    T.Horizon = [h];
    T.NumSamples = [summaryTrain.NumSamples];
    T.UniqueSites = [summaryTrain.UniqueSites];
    T.YearRange = sprintf('%d-%d',[summaryTrain.StartDate],[summaryTrain.EndDate]);
    T.avgDurationYears = [summaryTrain.DurationYears];
    T.MeanTarget = [summaryTrain.MeanTarget];
    T.StdTarget = [summaryTrain.StdTarget];

    summaryAll = [summaryAll; T];
end

for h = horizons
    fprintf('\n===== Horizon: +%d Years =====\n', h);

    testFile  = strjoin(["TestDatawERA5_", num2str(h), "Years_V3.csv"],'');

    % Read tables
    testData = readtable(testFile);

    % Summary for Test
    summaryTest = getSummary(testData, targetVar);

    % Combine and store
    T = table;
    T.Dataset = {'Test'};
    T.Horizon = [h];
    T.NumSamples = [summaryTest.NumSamples];
    T.UniqueSites = [summaryTest.UniqueSites];
    T.YearRange = sprintf('%d-%d',[summaryTest.StartDate],[summaryTest.EndDate]);
    T.avgDurationYears = [summaryTest.DurationYears];
    T.MeanTarget = [summaryTest.MeanTarget];
    T.StdTarget = [summaryTest.StdTarget];

    summaryAll = [summaryAll; T];
end
disp(summaryAll);

% features min-max ranges
clc
vrblz = trainData.Properties.VariableNames;
for ii = 4:82
    fprintf('%s, %d-%d\n',vrblz{ii},min([trainData{:,ii};testData{:,ii}]),max([trainData{:,ii};testData{:,ii}]))
end


%% Functions
function summary = getSummary(data, targetVar)
    summary.NumSamples = height(data);
    summary.UniqueSites = size(unique([data.Latitude data.Longitude], 'rows'), 1);
    summary.StartDate = min(data.Year);
    summary.EndDate = max(data.Year);

    if ismember(targetVar, data.Properties.VariableNames)
        [~,~,ic] = unique([data.Latitude data.Longitude], 'rows');
        for ii = 1:length(unique(ic))
            meanSite(ii) = mean(data{ic==ii,(targetVar)}, 'omitnan');
            stdSite(ii) = std(data{ic==ii,(targetVar)}, 'omitnan');
            durSite(ii) = max(data{ic==ii,'Year'})-min(data{ic==ii,'Year'});
        end
        summary.MeanTarget = mean(meanSite);
        summary.StdTarget = mean(stdSite);
        summary.DurationYears = mean(durSite);
    else
        summary.MeanTarget = NaN;
        summary.StdTarget = NaN;
        summary.DurationYears = NaN;        
    end
end

function T_final = innerJoinMultipleTables(tables)
    % Check input validity
    if isempty(tables)
        error('No tables provided!');
    end
    if ~iscell(tables)
        error('Input "tables" must be a cell array of tables.');
    end

    % Initialize with the first table
    T_final = tables{1};

    % Iteratively outer join with the remaining tables
    for k = 2:length(tables)
        T_final = innerjoin(T_final, tables{k}, "Keys",["Latitude", "Longitude","Year"]);
        T_final.Properties.VariableNames = strrep(T_final.Properties.VariableNames,'_T_final','');

        tmp = T_final.Properties.VariableNames(contains(T_final.Properties.VariableNames,'right')); %remove redundant columns of right table
        for iiii = 1:length(tmp)
            T_final.(tmp{iiii}) = [];
        end
    end
end

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

% Mann-Kendall Trend Strength
function tau = mann_kendall_trend(series)
    if length(series) < 3
        tau = NaN;
        return;
    end
    [tau, ~] = corr((1:length(series))', series, 'Type', 'Kendall', 'Rows', 'complete');
end

% Seasonality Strength (Fourier Analysis) with Error Handling
function strength = seasonality_strength_safe(series)
    if length(series) < 4  % Ensure enough data points
        strength = NaN;
        return;
    end
    freq_domain = abs(fft(series));
    if length(freq_domain) < 4  % Avoid indexing errors
        strength = NaN;
    else
        strength = sum(freq_domain(2:min(4, end))) / sum(freq_domain);
    end
end

% Lempel-Ziv Complexity
function lzc = lempel_ziv_complexity(series)
    if length(series) < 3
        lzc = NaN;
        return;
    end
    seq = normalize(series, 'range', [0, 255]);
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