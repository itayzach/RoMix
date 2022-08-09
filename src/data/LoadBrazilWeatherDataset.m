function sDataset = LoadBrazilWeatherDataset(sDatasetParams, N, n, nLabeled)
    T = readtable(fullfile('data','Brazilian_Weather_Stations-Temperature_1961-1990.xlsx'));
    
    % Remove entries with NaN
    T(isnan(T.Ann),:) = [];
    
    nPoints = height(T);
    latStr = cell2mat(T.Latitude);
    lonStr = cell2mat(T.Longitude);
    % Add seconds to lat/lon to have a valid str2angle input format:
    %     '07°38''S' --> '07°38''00"S'
    latStrFormatted = [latStr(:,1:5), repmat('''00"',nPoints,1), latStr(:,7)];
    latitude = str2angle(latStrFormatted);
    lonStrFormatted = [lonStr(:,1:5), repmat('''00"',nPoints,1), lonStr(:,7)];
    longitude = str2angle(lonStrFormatted);
    data = [longitude, latitude];
    
    trainTestRatio = n/N;
    N = length(data);
    n = min(round(trainTestRatio*N),N);
    N = length(data);
    rperm = randperm(N);
    dataRearranged = data(rperm,:);
    
    sDataset.sData.x = dataRearranged(1:n,:);
    sDataset.sData.xt = dataRearranged;
    
    nMonths = numel(sDatasetParams.xTickNames);
    mSignals = zeros(N, nMonths);
    for monthId = 1:nMonths
        currMonth = sDatasetParams.xTickNames{monthId};
        mSignals(:,monthId) = T.(currMonth)(rperm);
    end
    
    rpermUnlabeled = randperm(n);
    rpermUnlabeled = rpermUnlabeled(1:n-nLabeled);
    mSignalsMasked = mSignals(1:n,:);
    mSignalsMasked(rpermUnlabeled,:) = 0;

    sDataset.sData.y = mSignals(1:n,:);
    sDataset.sData.ymasked = mSignalsMasked;
    sDataset.sData.yt = mSignals;
end
