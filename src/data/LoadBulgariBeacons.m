function sDataset = LoadBulgariBeacons(sPlotParams, NLat, NLon, N, n, nLabeled)
assert(nLabeled <= n && n <= N);
%% Load DEM
fnameDem = fullfile('data', 'bulgari', 'DEM_Bulgari', 'N42E024_E24_N42.dt1');
[Z,R] = dted(fnameDem);
sFullMap = BuildFullMap(Z,R); %close(1);
nLon = size(sFullMap.lon,1);
nLat = size(sFullMap.lat,1);

%% Create grid
b_zoom = false;
if b_zoom
    %Build grid:
    latStartInd = 701;
    latEndInd = 800;

    lonStartInd = 451;
    lonEndInd = 550;
    dlat=(latEndInd-latStartInd+1)/NLat; 
    dlon=(lonEndInd-lonStartInd+1)/NLon;
    assert(dlat == round(dlat) && dlon == round(dlon), 'number of points must be integer')
    grid_lat_ind=latStartInd:dlat:latEndInd;
    grid_lon_ind=lonStartInd:dlon:lonEndInd;

    dn_lat = (nLat-2)/NLat; %jump in latitude every dn_lat points
    dn_lon = (nLon-2)/NLon; %jump in longitude every dn_lon points
    
    if sPlotParams.b_globalPlotEnable
        figure; ShowDEM(sFullMap.Z, sFullMap.R); hold on; 
        patch(sFullMap.lon([grid_lon_ind(1)*ones(1,2),grid_lon_ind(end)*ones(1,2)]), ...
              sFullMap.lat([grid_lat_ind(1),grid_lat_ind(end)*ones(1,2),grid_lat_ind(1)]),...
              'b:', 'FaceColor','none','EdgeColor','blue');

        %grid=BuildGrid(grid_lat_ind, grid_lon_ind, sFullMap);   
        %figure; ShowDEM(grid.Z, grid.R);
    end
else
    dn_lat = (nLat-2)/NLat; %jump in latitude every dn_lat points
    dn_lon = (nLon-2)/NLon; %jump in longitude every dn_lon points
    grid_lat_ind = dn_lat/2+1:dn_lat:nLat-dn_lat/2;
    grid_lon_ind = dn_lon/2+1:dn_lon:nLon-dn_lon/2;
    assert((nLat-2)/sqrt(N) == floor((nLat-2)/sqrt(N)) && (nLon-2)/sqrt(N) == floor((nLon-2)/sqrt(N)))
end

%grid=BuildGrid(grid_lat_ind, grid_lon_ind, sFullMap);

nGridLon = size(grid_lon_ind,2);
nGridLat = size(grid_lat_ind,2);
vGridLon = sFullMap.lon(grid_lon_ind);
vGridLat = sFullMap.lat(grid_lat_ind);
[mGridLon, mGridLat] = meshgrid(vGridLon, vGridLat);
mGridZ = double(sFullMap.Z(grid_lat_ind, grid_lon_ind));
if b_zoom
    vGridR = [(nGridLat-1)/(vGridLat(end)-vGridLat(1)), vGridLat(end), vGridLon(1)]; %redefine the number of cells
    if sPlotParams.b_globalPlotEnable
        figure; ShowDEM(mGridZ,vGridR);
    end
else
    vGridR = [nGridLat-1, sFullMap.R(2), sFullMap.R(3)];
end

Re=6371000; %[m] earth radius.
fprintf(1,'Using grid of %d x %d points (resolution %.1f[Km] x %.1f[Km])\n',...
    length(grid_lat_ind),length(grid_lon_ind), ...
    (2*pi*Re/360/1e3)*diff(vGridLat(1:2)),(2*pi*Re/360/1e3)*diff(vGridLon(1:2)));
fprintf(1,'Total tested area of %.1f[Km] x %.1f[Km]\n', (2*pi*Re/360/1e3)*(vGridLat(end)-vGridLat(1)),(2*pi*Re/360/1e3)*(vGridLon(end)-vGridLon(1)));

assert(numel(mGridZ(:)) == N)

%% Measurments
NoiseFloor_dBm = -154; 
noise_level_mW = 10^(NoiseFloor_dBm/10);
[tx, sense, y, y_M] = MyBuildMeasurements(sPlotParams, nLabeled, vGridLat, vGridLon, noise_level_mW);

%% Plot
if sPlotParams.b_globalPlotEnable
    for i=1:2
        if i == 1
            figName = 'BulgariDEM_noGrid';
        elseif i == 2
            figName = 'BulgariDEM_grid';
        end
        fig1 = figure('Name', figName);
        ShowDEM(sFullMap.Z, sFullMap.R); hold on;
        %ShowDEM(mGridZ, vGridR); hold on;
        plot(sense.lon, sense.lat, 'yo','MarkerFaceColor','y','MarkerSize',3);
        %text(sense.lon, sense.lat, num2str([1:Nsensors]'));
        plot(tx.lon, tx.lat, 'bd','MarkerFaceColor','b','MarkerSize',3); %test sites
        if i == 2
            plot(mGridLon, mGridLat,'k.','MarkerSize',2);
        end
        if b_zoom
            patch(sFullMap.lon([grid_lon_ind(1)*ones(1,2),grid_lon_ind(end)*ones(1,2)]),...
                  sFullMap.lat([grid_lat_ind(1),grid_lat_ind(end)*ones(1,2),grid_lat_ind(1)]),...
                  'b:', 'FaceColor','none','EdgeColor','blue');
        end
        set(gca,'FontSize', 14);
        x0 = 10; y0 = 50; height = 400; width = 500;
        set(gcf,'Position', [x0 y0 width height])
        SaveFigure(sPlotParams, fig1, figName, {'epsc', 'png'});
    end
    
    % Just for visualization of missing points in the DEM
    rperm2 = randperm(numel(sFullMap.Z(:)));
    mFullGridZmissing = sFullMap.Z;
    mFullGridZmissing(rperm2(1:100000)) = nan;
    vGridR = sFullMap.R;
    
    figName = 'BulgariDEM_missing';
    fig1 = figure('Name', figName);
    ShowDEM(mFullGridZmissing, vGridR); hold on;
    plot(sense.lon, sense.lat, 'yo','MarkerFaceColor','y','MarkerSize',3);
    plot(tx.lon, tx.lat, 'bd','MarkerFaceColor','b','MarkerSize',3); %test sites
    set(gca,'FontSize', 14);
    x0 = 10; y0 = 50; height = 400; width = 500;
    set(gcf,'Position', [x0 y0 width height])
    SaveFigure(sPlotParams, fig1, figName, {'epsc', 'png'});
end


%% Visibility (LoS)
if b_zoom
    visPath = fullfile('data', 'bulgari', ['tVis_zoom_', 'lat', num2str(nGridLat) '_lon', num2str(nGridLon) '.mat']);
else
    visPath = fullfile('data', 'bulgari', ['tVis_', 'lat', num2str(nGridLat) '_lon', num2str(nGridLon) '.mat']);
end
if ~isfile(visPath)
   assert(all(~isnan(mGridZ(:))))
   tVis = MyCalcVisibility(mGridZ, vGridR, vGridLat, vGridLon);
   save(visPath,'tVis')
else
   load(visPath,'tVis')
end
A_LoS = zeros(size(tVis,3));
N = nGridLon*nGridLat; % number of grid points
for i = 1:N
   A_LoS(:,i) = reshape(tVis(:,:,i), N, []);
end
A_LoS = max(A_LoS,A_LoS'); % assure symmetry (assumes LoS if one way exists)
A2_LoS = A_LoS*A_LoS; % 2-hop neighborhoods

%% Free Space Attenuation
wgs84 = wgs84Ellipsoid('meter');
[X,Y,Z] = geodetic2ecef(wgs84,mGridLat(:),mGridLon(:),mGridZ(:));

% Due to the grid sampling, we get distance = 0 on the diagonal.
% To fix this, we set the distance of node to itself according to farthest point in the cell: 
sub_grid_distance = 0.5*Re*diff(sFullMap.lat(1:2,1))*(pi/180)*dn_lat;

dist = squareform(pdist([X(:), Y(:), Z(:)])) + sub_grid_distance*eye(N);
fRF = 150e6; c = 3e8; lambda = c/fRF;
FSPL = (4*pi*dist/lambda).^2;

%% Build sDataset
latent = [X, Y, Z, FSPL(:,tx.ind), A_LoS(:,tx.ind), A2_LoS(tx.ind,:).'];
latent = (latent - min(latent))./(max(latent) - min(latent));
assert(size(latent,2) == sPlotParams.dim);

%% Save dataset and grid for plots
sDataset = BuildDataset(N, n, latent, sense, tx, y, y_M, mGridLon, mGridLat, vGridLon, vGridLat);
rng(3);
% grid.lat = sDataset.sData.vGridLat;
% grid.lon = sDataset.sData.vGridLon;
% mGridShuffledInt = sDataset.sData.mGrid(sDataset.sData.vShuffleMap,:);
% mGridShuffledRec = mGridShuffledInt(1:n,:);
% 
% ytOrig = zeros(N,1);
% ytOrig(sDataset.sData.vShuffleMap) = sDataset.sData.yt;
% vUnshuffledSigIntRef = ytOrig;
% PlotGraphSignals(sPlotParams, [], 'GroundTruthInterp', ...
%    {mGridShuffledRec}, {sDataset.sData.y}, {'$s$'}, {1:n}, ...
%    [], [], [], [], grid, sDataset.sData.tx, sDataset.sData.sense);
% PlotGraphSignals(sPlotParams, [], 'GroundTruthExtrap', ...
%    {sDataset.sData.mGrid}, {vUnshuffledSigIntRef}, {'$\tilde{s}$'}, {1:N}, ...
%    [], [], [], [], grid, sDataset.sData.tx, sDataset.sData.sense);


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
    fig2 = figure('Name', 'Signal-pro'); 
    colormap(jet);
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
    set(gca,'FontSize', 14);
    x0 = 10; y0 = 50; height = 400; width = 500;
    set(gcf,'Position', [x0 y0 width height])
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
%figure;(1); clf; %h=ShowDEM(Z,R);
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

fig = figure;
h=ShowDEM(fmap.Z,fmap.R);
fmap.Nlat=size(Z,1); %# of latitudes entries
fmap.Nlon=size(Z,2); %# of longitudes entries
fmap.lon=h.XData(1,:)'; %all lon points in the figure
fmap.lat=h.YData(:,1); %all lat points in the figure
close(fig);

end

function h=ShowDEM(Z,R) %show elevation map
h=geoshow(Z,R,'DisplayType','texturemap','CData',Z);
zlimits = [min(Z(:)) max(Z(:))];
demcmap(zlimits); %Colormaps appropriate to terrain elevation data
c=colorbar;  
c.Label.String='Height [m]'; 
c.Label.Interpreter='latex'; 
c.Label.FontSize = 14;
c.TickLabelInterpreter = 'latex';
axis tight;
xlabel('Lon. [deg]','FontSize', 14);
ylabel('Lat. [deg]','FontSize', 14);
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


function sDataset = BuildDataset(N, n, latent, sense, tx, y, y_M, mGridLon, mGridLat, vGridLon, vGridLat)
% Arrange data to 
% 1     ->  ell : sensors (measurements)
% ell+1 ->  n   : for interp
% n+1   ->  N   : for extrap
ell = numel(sense.ind');
vIndNoSense = setdiff(1:N,sense.ind');
rperm = randperm(N-ell);
vIndNoSenseShuffled = vIndNoSense(rperm);

vShuffleMap = [sense.ind, vIndNoSenseShuffled]';

dataRearranged = latent(vShuffleMap,:);
yRearranged = y(vShuffleMap);
assert(isequal(yRearranged(1:ell), y_M))
sDataset.sData.x = dataRearranged(1:n,:);
sDataset.sData.xt = dataRearranged;
sDataset.sData.y = 10*log10(abs(yRearranged(1:n))); % from Signal-Pro 
%sDataset.sData.y = yRearranged(1:n); % from Signal-Pro 
sDataset.sData.ymasked = zeros(n,1);
sDataset.sData.ymasked(1:ell) = 10*log10(abs(y_M)); % vLabeledInd = sense.ind';
%sDataset.sData.ymasked(1:ell) = y_M; % vLabeledInd = sense.ind';
sDataset.sData.yt = 10*log10(abs(yRearranged)); % from Signal-Pro 
%sDataset.sData.yt = yRearranged; %from Signal-Pro 

% Used only for plots
sDataset.sData.mGrid = [mGridLon(:), mGridLat(:)];
sDataset.sData.vGridLon = vGridLon;
sDataset.sData.vGridLat = vGridLat;
sDataset.sData.sense = sense;
sDataset.sData.tx = tx;
sDataset.sData.vShuffleMap = vShuffleMap;

VerifyShuffle(sDataset, latent, N, n, ell, sense, y, y_M, vShuffleMap);
end


function VerifyShuffle(sDataset, data, N, n, ell, sense, y, y_M, vShuffleMap)
dim = size(data,2);
xtOrig = zeros(N,dim);
xtOrig(vShuffleMap,:) = sDataset.sData.xt;
ytOrig = zeros(N,1);
ytOrig(vShuffleMap) = sDataset.sData.yt;
xOrig = xtOrig(1:n,:);
yOrig = ytOrig((1:n)');
ymaskedOrig = zeros(n,1);
ymaskedOrig((1:ell)') = ytOrig(sense.ind');

assert(isequal(xtOrig, data))
assert(isequal(xOrig, data(1:n,:)))
%assert(isequal(yOrig, y(1:n)))
assert(isequal(yOrig, 10*log10(abs(y(1:n)))))
%assert(isequal(ymaskedOrig(1:ell), y_M))
assert(isequal(ymaskedOrig(1:ell), 10*log10(abs(y_M))))
assert(isequal(ymaskedOrig(ell+1:end), zeros(n-ell,1)))
%assert(isequal(ytOrig, y)) 
assert(isequal(ytOrig, 10*log10(abs(y))))
end