function sDataset = LoadBulgariBeacons(sPlotParams, NLat, NLon, N, n, nLabeled)
%% Load DEM
fnameDem = fullfile('data', 'bulgari', 'DEM_Bulgari', 'N42E024_E24_N42.dt1');
[Z,R] = dted(fnameDem);
sFullMap = BuildFullMap(Z,R); %close(1);
LON = sFullMap.lon;
LAT = sFullMap.lat;

nLon = size(LON,1);
nLat = size(LAT,1);

%% Create grid
assert((nLat-2)/sqrt(N) == floor((nLat-2)/sqrt(N)) && (nLon-2)/sqrt(N) == floor((nLon-2)/sqrt(N)))
b_zoom = false;
if b_zoom
    % Not working so well...
    dlat = 100/NLat; %jump in latitude every dn_lat points
    dlon = 200/NLon; %jump in longitude every dn_lon points
    grid_lat_ind=100:dlat:200-1;
    grid_lon_ind=1:dlon:200;
else
    dn_lat = (nLat-2)/NLat; %jump in latitude every dn_lat points
    dn_lon = (nLon-2)/NLon; %jump in longitude every dn_lon points
    grid_lat_ind = dn_lat/2+1:dn_lat:nLat-dn_lat/2;
    grid_lon_ind = dn_lon/2+1:dn_lon:nLon-dn_lon/2;
end

%grid=BuildGrid(grid_lat_ind, grid_lon_ind, sFullMap);

nGridLon = size(grid_lon_ind,2);
nGridLat = size(grid_lat_ind,2);
vGridLon = LON(grid_lon_ind);
vGridLat = LAT(grid_lat_ind);
[mGridLon, mGridLat] = meshgrid(vGridLon, vGridLat);
mGridZ = double(sFullMap.Z(grid_lat_ind, grid_lon_ind));
if b_zoom
    vGridR = [(nGridLat-1)/(vGridLat(end)-vGridLat(1)), vGridLat(end), vGridLon(1)]; %redefine the number of cells
    figure; ShowDEM(mGridZ(:),vGridR);
else
    vGridR = [nGridLat-1, sFullMap.R(2), sFullMap.R(3)];
end

Re=6371000; %[m] earth radius.
fprintf(1,'Using grid of %d x %d points (resolution %.1f[Km] x %.1f[Km])\n',...
    length(grid_lat_ind),length(grid_lon_ind), ...
    (2*pi*Re/360/1e3)*diff(vGridLat(1:2)),(2*pi*Re/360/1e3)*diff(vGridLon(1:2)));
fprintf(1,'Total tested area of %.1f[Km] x %.1f[Km]\n', (2*pi*Re/360/1e3)*(vGridLat(end)-vGridLat(1)),(2*pi*Re/360/1e3)*(vGridLon(end)-vGridLon(1)));

assert(numel(mGridZ(:)) == N)
%% Visibility (LoS)
visPath = fullfile('data', 'bulgari', ['tVis_', 'lat', num2str(nGridLat) '_lon', num2str(nGridLon) '.mat']);
if ~isfile(visPath)
   assert(all(~isnan(mGridZ(:))))
   tVis = MyCalcVisibility(mGridZ, vGridR, vGridLat, vGridLon);
   save(visPath,'tVis')
else
   load(visPath,'tVis')
end
A_LoS_2m = zeros(size(tVis,3));
N = nGridLon*nGridLat; % number of grid points
for i = 1:N
   A_LoS_2m(:,i) = reshape(tVis(:,:,i), N, []);
end
A_LoS_2m = max(A_LoS_2m,A_LoS_2m'); % assure symmetry (assumes LoS if one way exists)
A_LoS_2m = A_LoS_2m-diag(diag(A_LoS_2m)); % remove self loops
A_LoS_2m_sq = A_LoS_2m*A_LoS_2m; % 2-hop neighborhoods
%A_LoS_2m_qub = A_LoS_2m*A_LoS_2m_sq; % 2-hop neighborhoods

%% Free Space Attenuation
attenuationAlpha = 2;
WITH_SELF_LOOPS = false;
Atten2 = MyCalcAttenuation(mGridLat, mGridLon, mGridZ, attenuationAlpha, WITH_SELF_LOOPS);
A_gain_mW = 10.^(-Atten2/10);

%Fix synthesis model to describe small attenuation (large gain) from a node to iself.
Re=6371000; %[m] earth radius.
sub_grid_distance = 0.5*Re*diff(LAT(1:2,1))*(pi/180)*dn_lat;
sub_grid_gain = 1/(sub_grid_distance^2);
A_gain_mW = A_gain_mW + sub_grid_gain*eye(N); %gain only mat
A_gain_dB = 10*log10(A_gain_mW);

%% Measurments
NoiseFloor_dBm = -154; 
noise_level_mW = 10^(NoiseFloor_dBm/10);
[tx, sense, y, y_M] = MyBuildMeasurements(sPlotParams, nLabeled, vGridLat, vGridLon, noise_level_mW);

%% Build sDataset
assert(N == n);
assert(nLabeled <= n);
wgs84 = wgs84Ellipsoid('kilometer');
[X,Y,Z] = geodetic2ecef(wgs84,mGridLat(:),mGridLon(:),mGridZ(:));

A_LoS_2m_plusMinus = (-1+2*A_LoS_2m);
b_normalizeNormal = false;
b_normalizeUniform = false;

% data = [X, Y, Z]; b_normalizeUniform = true;
% data = [X, Y, Z, A_gain_dB(:,tx.ind)];  b_normalizeUniform = true;%
% data = [X, Y, Z, A_LoS_2m(:,tx.ind)];  b_normalizeUniform = true;%
% data = [X, Y, Z, A_gain_dB(:,tx.ind), A_LoS_2m(:,tx.ind)];  b_normalizeUniform = true;% 92.34% omega = 1.5; gmm = 10; gamma = 0.05; nLabeled = 50;
% data = [Z, A_gain_dB(:,tx.ind), A_LoS_2m(:,tx.ind), A_LoS_2m_sq(:,tx.ind)];  b_normalizeUniform = true; % 92.56%  omega = 1.5; gmm = 10; gamma = 0.05; nLabeled = 50;
data = [X, Y, Z, A_gain_dB(:,tx.ind), A_LoS_2m(:,tx.ind), A_LoS_2m_sq(:,tx.ind)];  b_normalizeUniform = true; % 92.45% omega = 1.5; gmm = 10; gamma = 0.05; nLabeled = 50;


% Normalize
if b_normalizeUniform
    data = (data - min(data))./(max(data) - min(data));
end
if b_normalizeNormal    
    data = (data - mean(data)) ./ std(data);
end

%% Plot
if sPlotParams.b_globalPlotEnable
    fig1 = figure(1); clf;
    ShowDEM(sFullMap.Z, sFullMap.R); hold on;
    plot(sense.lon, sense.lat, 'yo','MarkerFaceColor','y','MarkerSize',3);
    %text(sense.lon, sense.lat, num2str([1:Nsensors]'));
    plot(tx.lon, tx.lat, 'bd','MarkerFaceColor','b','MarkerSize',3); %test sites
    plot(mGridLon, mGridLat,'k.','MarkerSize',2);
    assert(size(data,2) == 9);
    SaveFigure(sPlotParams, fig1, 'BulgariDEM', {'epsc', 'png'});
end

%% Save dataset and grid for plots
sDataset.sData.x = data;
sDataset.sData.xt = data;
sDataset.sData.y = 10*log10(abs(y)); % from Signal-Pro 
sDataset.sData.ymasked = zeros(n,1);
sDataset.sData.ymasked(sense.ind') = 10*log10(abs(y_M)); % vLabeledInd = sense.ind';
sDataset.sData.yt = 10*log10(abs(y)); % from Signal-Pro 

% Used only for plots
sDataset.sData.mGrid = [mGridLon(:), mGridLat(:)];
sDataset.sData.vGridLon = vGridLon;
sDataset.sData.vGridLat = vGridLat;
sDataset.sData.sense = sense;
sDataset.sData.tx = tx;


end
%% Local functions
% =============================================================================================
%Collect measurements from the sigpro (or from ideal model)
function [tx, sense, y, y_M]=MyBuildMeasurements(sPlotParams, Nsensors, vGridLat, vGridLon, noise_level_mW)

grid.lat = vGridLat;
grid.lon = vGridLon;
grid.Nlat = numel(vGridLat);
grid.Nlon = numel(vGridLon);

%Load Edi's simulation vectors:
fprintf(1,'Reading Signal-Pro results... ');
[sensors_lat, sensors_lon, sensors_alt]= GetSensorsPos();
simul_path=fullfile('data', 'bulgari', 'SensorsSimulation');
for siteInd=1:5
    load(fullfile(simul_path,sprintf('site%d.mat',siteInd)));
    Site(siteInd).lat=site.lat;
    Site(siteInd).lon=site.lon;
    Site(siteInd).Atten=-site.PL;    
    Site(siteInd).LoS=site.LoS;       
    Site(siteInd).tx_lat=sensors_lat(siteInd);
    Site(siteInd).tx_lon=sensors_lon(siteInd);
    
    %receiver/transmistter but sampled on the grid:
    %[rt(active_sites).Atten, rt(active_sites).LoS, rt(active_sites).grid_lat_ind, rt(active_sites).grid_lon_ind]=site2grid(Site(active_sites), grid); %rt for receive/transmit
    [rt(siteInd).Atten, rt(siteInd).LoS, rt(siteInd).grid_lat_ind, rt(siteInd).grid_lon_ind]=MyConvSiteToGrid(Site(siteInd), grid);
    
    rt(siteInd).lat=sensors_lat(siteInd);
    rt(siteInd).lon=sensors_lon(siteInd);
    rt(siteInd).Power_dBm=0;
end
fprintf(1,'Done\n');

%Set transmitters
%----------------
active_sites=[4:5]; %the active sites to work with
rt = rt(active_sites);
%Define the transmitter by the selected sites: 
tx.lat=[rt.lat]; 
tx.lon=[rt.lon];
Ns=length(rt);
tx.alt=2*ones(1,Ns); %[m] Above map
tx.Power_dBm=0*ones(1,Ns);
[tx.grid_lat, tx.grid_lon, tx.grid_lat_ind, tx.grid_lon_ind, tx.ind]=findClosest(tx.lat, tx.lon, grid); %closest position on the grid

%Set sensors
%-----------
SENSORS_ON_GRID=true;
sense=SetSensors(Nsensors, grid, tx, SENSORS_ON_GRID);
%shift one of the sensors to make the problem less trivial (it was too close
%to a source)
% sense.lat_ind(20)=sense.lat_ind(20)+3; sense.lat(20)=grid.lat(sense.lat_ind(20));
% sense.lon_ind(20)=sense.lon_ind(20)-5; sense.lon(20)=grid.lon(sense.lon_ind(20));
% sense.ind(20)=sub2ind([grid.Nlat, grid.Nlon], sense.lat_ind(20), sense.lon_ind(20));
%warning('In Tsvi''s code the following 3 lines are without comment')
%sense.ind(20)=830; [sense.lat_ind(20),sense.lon_ind(20)]=ind2sub([grid.Nlat, grid.Nlon],sense.ind(20));
%sense.lat(20)=grid.lat(sense.lat_ind(20));
%sense.lon(20)=grid.lon(sense.lon_ind(20));



%Build measurements:
%-------------------
y=0;
if 0 %params.IDEAL_MODEL    
    for(ns=1:Ns) %go over all sources
        y_source_mw=Amodel(:,tx.ind(ns))*ns; %Full grid. In dBm
        y=y+y_source_mw;
    end
else
    for(ns=1:Ns) %go over all sources
        y_source_dBm=rt(ns).Atten(:); %Full grid. In dBm        
        
        %add up powers from different sources
        y=y+10.^(y_source_dBm/10); %full grid vector. (Reception power). In mw
    end
    %set some saturation due to noise:
    y=max(y,noise_level_mW);
end

%nGridPoints = grid.Nlat*grid.Nlon;
%I=eye(nGridPoints); S_M=I(sense.ind,:);
%y_M=S_M*y; % sampling matrix. Same as y_M = y(sense.ind)
y_M = y(sense.ind);

%Show the signal pro map
if sPlotParams.b_globalPlotEnable
    fig2 = figure(2); clf; colormap(hot);   
    %calc total received power due to SignalPro:
    y_sigPro_mw=0;
    for(source=1:Ns)
        site_ind=active_sites(source);
        y_sigPro_mw=y_sigPro_mw+10.^((Site(site_ind).Atten(:)+tx.Power_dBm(source))/10);
    end  
    y_sigPro_mw=max(y_sigPro_mw,noise_level_mW);
    PlotPowers(10*log10(y_sigPro_mw), Site(1)); %title('Signal Pro');
    hold on;
    plot(tx.lon, tx.lat, 'db','MarkerFaceColor','b','MarkerSize',2.5);
    plot(sense.lon, sense.lat, 'yo','MarkerSize',3);
    SaveFigure(sPlotParams, fig2, 'BulgariSignalPro', {'epsc', 'png'});
end
end

function [grid_lat, grid_lon, grid_lat_ind, grid_lon_ind, grid_ind]=findClosest(lat, lon, grid)
[~,pos_lat]=min(abs(lat-grid.lat));
[~,pos_lon]=min(abs(lon-grid.lon));
grid_lat=grid.lat(pos_lat);
grid_lon=grid.lon(pos_lon);
grid_lat_ind=pos_lat;
grid_lon_ind=pos_lon;
grid_ind=sub2ind([grid.Nlat, grid.Nlon],grid_lat_ind,grid_lon_ind);
end

function sense=SetSensors(Nsensors, grid, tx, SENSORS_ON_GRID)

nGridPoints = grid.Nlat*grid.Nlon;

%Define position of sensors:
rng(3); %for reproducability
%warning('I comment rng')

Re=6371000; %[m] earth radius.
MIN_ALLOWED_KM = 5.5;
MIN_ALLOWED_DEG = MIN_ALLOWED_KM/((2*pi*Re/360/1e3));
mDistToSense = zeros(Nsensors,numel(tx.alt));
iter = 0;
fprintf('locating sensors..')
while sum(mDistToSense(:) < MIN_ALLOWED_DEG)
    if SENSORS_ON_GRID
        sense.ind=randperm(nGridPoints,Nsensors); %random columns
        %[sense.lat_ind,sense.lon_ind]=ind2sub([grid.Nlat, grid.Nlon],sense.ind);
        [sense.lat_ind,sense.lon_ind]=ind2sub([grid.Nlat, grid.Nlon],sense.ind);
        sense.lat=grid.lat(sense.lat_ind);
        sense.lon=grid.lon(sense.lon_ind);
    else
        sense.lat=grid.lat(1)+rand(Nsensors,1)*(grid.lat(end)-grid.lat(1));
        sense.lon=grid.lon(1)+rand(Nsensors,1)*(grid.lon(end)-grid.lon(1));
        sense.lat_ind=zeros(Nsensors,1);
        sense.lon_ind=zeros(Nsensors,1);
        for(n=1:Nsensors)
            [~,sense.lat_ind(n)] = min(abs(grid.lat-sense.lat(n)));
            [~,sense.lon_ind(n)] = min(abs(grid.lon-sense.lon(n)));
        end
        sense.ind=sub2ind([grid.Nlat, grid.Nlon],sense.lat_ind,sense.lon_ind);
    end
    
    for txInd = 1:numel(tx.alt)
        mDistToSense(:,txInd) = vecnorm([tx.lat(txInd), tx.lon(txInd)] - [sense.lat, sense.lon],2,2);    
    end
    iter = iter + 1;
    fprintf('.');
end
fprintf(' took %d iterations.\n', iter)

end

function W = MyCalcAttenuation(Lat, Lon, Z, attenuation_alpha, WITH_SELF_LOOPS)

if nargin<5 || isempty(WITH_SELF_LOOPS)
    WITH_SELF_LOOPS=false; %by default, assume inifinte attenuation from a point to itself, to make this equivalent to no self loops
end

if nargin<4 || isempty(attenuation_alpha)
    attenuation_alpha=2; %assume attenuation proportional to square of the distance (free space model) by default.
end

%Transform horizontal position to meters:
ref_lat=Lat(1); ref_lon=Lon(1); %use this as a reference point
Re=6371000; %[m] earth radius.
%This is not exact but hopefully good enought. Go from degree difference,
%to meters difference:
% X,Y are not actually X,Y! They are okay only for distance calculation because some of the entires are 0!!!

X=(Lat(:)-ref_lat)*(pi/180)*Re;
Y=(Lon(:)-ref_lon)*(pi/180)*Re;

%wgs84 = wgs84Ellipsoid('kilometer');
%[X,Y,Z] = geodetic2ecef(wgs84,Lat(:),Lon(:),Z(:));


pd=pdist([X(:),Y(:),Z(:)], 'euclidean'); %all pairs of euclidean distances
atten=10*log10(pd.^attenuation_alpha); %attenuation proportional to square of the distance in free space
W = squareform(atten);
if WITH_SELF_LOOPS==false
    W = W + 300*eye(size(W)); %plug large number along the diagonal to set large attenuation from a node to itself
end
end

function tVis = MyCalcVisibility(Z, R, vGridLat, vGridLon)

rxAlt = 2; %[m]
txAlt = 2; %[m]

nLat = numel(vGridLat);
nLon = numel(vGridLon);

tVis = zeros(nLat,nLon,nLat*nLon); %malloc
fprintf(1,'calculating LoS map for each measuerement point on the grid:      ');
for n = 1:nLat
    for m = 1:nLon
        ind=sub2ind([nLat,nLon],n,m);
        tVis(:,:,ind)=viewshed(Z, R, vGridLat(n), vGridLon(m), rxAlt, txAlt);
        fprintf(1,'\b\b\b\b\b');
        fprintf(1,'%4.1f%%',100*((n-1)*nLon+m)/(nLat*nLon));
    end
end
fprintf(1,' Done.\n');
end

function [PathLoss, LoS, tx_grid_lat_ind, tx_grid_lon_ind]=MyConvSiteToGrid(site, grid)
lat_ind=zeros(grid.Nlat,1); %to store positions of Edi's site data along latitude which correspond best to the grid lat
for(n=1:grid.Nlat)
    [~,lat_ind(n)]=min(abs(grid.lat(n)-site.lat));
end
lon_ind=zeros(grid.Nlon,1);
for(n=1:grid.Nlon)
    [~,lon_ind(n)]=min(abs(grid.lon(n)-site.lon));
end
[LAT_IND,LON_IND]=meshgrid(lat_ind, lon_ind);
PathLoss=site.Atten(lat_ind,lon_ind);
[~,p]=hist(site.LoS,2); %the values are not around 0 and 1, so fixing this.
LoS=site.LoS(lat_ind,lon_ind);
LoS(LoS>mean(p))=1;
LoS(LoS<mean(p))=0;
%figure; subplot(211); imagesc(site.LoS); subplot(212); imagesc(LoS);

%Writing position of transmitter in grid terms:
[~,tx_grid_lat_ind]=min(abs(grid.lat-site.tx_lat));
[~,tx_grid_lon_ind]=min(abs(grid.lon-site.tx_lon));
end

function [sensors_lat, sensors_lon, sensors_alt]= GetSensorsPos()

data= [...
42	46	29.33	24	10	14.17	1533.549072
42	40	45.68	24	38	30.52	453.4824219
42	20	40.38	24	29	0.86	272.8933105
42	20	8.6	    24	46	56.25	201.5813904
42	35	50.14	24	24	37.27	1211.85437
42	23	28.35	24	10	13.54	793.7153931
42	48	57.93	24	49	17.54	1279.878784
42	11	23.65	24	30	57.7	195.7651367
42	11	1.72	24	54	28.85	153.1878815
42	8	44.2	24	7	56.6	404.8523865
42	42	58.54	24	55	26.96	2175.799072
42	45	7.01	24	24	42.21	2114.049805
42	21	56.04	24	5	39.89	329.4684448
42	37	30.41	24	16	9.34	1344.955078
42	24	5.71	24	9	44.35	768.0029297
42	32	6.29	24	42	15.19	468.918396
42	54	17.67	24	19	26.76	618.5600586
42	53	4.58	24	39	13.87	538.5880737
42	35	8.17	24	46	57.86	317.9203796
42	14	47.5	24	41	53.01	189.9951782];

sensors_lat=data(:,1)+data(:,2)/60+data(:,3)/3600;
sensors_lon=data(:,4)+data(:,5)/60+data(:,6)/3600;
sensors_alt=data(:,7);
end

function fmap = BuildFullMap(Z,R)
figure(1); clf; %h=ShowDEM(Z,R);
% if PRINT
%     fname='pics\Bulgari'; print(fname,'-dpdf'); print(fname,'-depsc'); print(fname,'-dtiff');
% end


%full map:
fmap.Z=Z; fmap.R=R;
%Fix problematic NaN values within Z
bad_ind=find(isnan(fmap.Z));
for(n=1:length(bad_ind))
    ind_start=max(1,bad_ind(n)-1);
    ind_stop=min(bad_ind(n)+1, numel(fmap.Z));
    fmap.Z(bad_ind(n))=mean(max(fmap.Z(ind_start:ind_stop)),min(fmap.Z(ind_start:ind_stop)));
end
if any(isnan(fmap.Z(:)))
    warning('There are still NaN within Z');
end

h=ShowDEM(fmap.Z,fmap.R);
    
fmap.Nlat=size(Z,1); %# of latitudes entries
fmap.Nlon=size(Z,2); %# of longitudes entries
fmap.lon=h.XData(1,:)'; %all lon points in the figure
fmap.lat=h.YData(:,1); %all lat points in the figure

end

function h=ShowDEM(Z,R) %show elevation map
h=geoshow(Z,R,'DisplayType','texturemap','CData',Z);
zlimits = [min(Z(:)) max(Z(:))];
demcmap(zlimits); %Colormaps appropriate to terrain elevation data
c=colorbar;  c.Label.String='Height [m]';
axis tight;
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
end

function grid=BuildGrid(grid_lat_ind, grid_lon_ind, fmap)
grid.lat_ind=grid_lat_ind; %indices with respect to the full map
grid.lon_ind=grid_lon_ind; %indices with respect to the full map
grid.lat=fmap.lat(grid_lat_ind); grid.Nlat=length(grid.lat_ind);
grid.lon=fmap.lon(grid_lon_ind); grid.Nlon=length(grid.lon_ind);
[gridLat,gridLon]=meshgrid(grid.lat, grid.lon); %grid on the map axis (lat lon)
grid.LAT=gridLat;
grid.LON=gridLon;
grid.Z=fmap.Z(grid_lat_ind, grid_lon_ind);
grid.R=[(grid.Nlat-1)/(grid.lat(end)-grid.lat(1)), grid.lat(end), grid.lon(1)]; %redefine the number of cells
end