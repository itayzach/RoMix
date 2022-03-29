function dist = CalcDistance(x1, x2, sDistanceParams)
if strcmp(sDistanceParams.distType, 'Euclidean')
    dist = pdist2(x1, x2);
elseif strcmp(sDistanceParams.distType, 'Haversine')
    earthRadius = sDistanceParams.earthRadius;
    lat1 = deg2rad(x1(:,1));
    lon1 = deg2rad(x1(:,2));
    lat2 = deg2rad(x2(:,1));
    lon2 = deg2rad(x2(:,2));
    
    [LON2,LON1] = meshgrid(lon2,lon1);
    [LAT2,LAT1] = meshgrid(lat2,lat1);
    DLAT = LAT2 - LAT1;
    DLON = LON2 - LON1;
    A = sin(DLAT/2).^2 + cos(LAT1).*cos(LAT2).*(sin(DLON/2).^2);
    dist = 2 * earthRadius * asin(sqrt(A));
  
else
    error('invalid distType')
end

end