% Data from: ftp://ftp.ncdc.noaa.gov/pub/data/gsod/
%% Restart
clc; clear; close all;
rng('default');

%% Read weather table
dataFile = fullfile('data', 'temperature', 'isd-history.csv');
opts = detectImportOptions(dataFile);
opts.SelectedVariableNames = {...
    'USAF', ...          Station number
    'STATIONNAME', ...   Station name
    'CTRY', ...          Country
    'LAT', ...           Latitude  [-90 deg, 90 deg]
    'LON', ...           Longitude [-180 deg, 180 deg]
    'BEGIN', ...         Begin date: 20010807 is 07/08/2001 (day/month/year)
    'END' ...            End date
    };
    %'ELEV_M_', ...       Elevation [meters]

opts.MissingRule = 'omitrow';
weatherTable = readtable(dataFile, opts);

%% Get only USA stations with data for 01/01/2001
desiredDate = 20100101;
missingIdentification = -999;
usaWeatherRows = (strcmpi(weatherTable.CTRY,'US')& ...
                  weatherTable.BEGIN <= desiredDate & ...
                  weatherTable.END >= desiredDate);
                  %weatherTable.ELEV_M_ ~=  missingIdentification);% & ...              
usaWeatherTable = weatherTable(usaWeatherRows,:);
% summary(usaWeatherTable);

%% Plot data
figure;
sz = 5;
scatter(weatherTable.LON, weatherTable.LAT, sz, 'filled');

GMModel = fitgmdist([usaWeatherTable.LON usaWeatherTable.LAT],1);

