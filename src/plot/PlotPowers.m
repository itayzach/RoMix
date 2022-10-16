function PlotPowers(power_vec, grid, tick_vals, SHOW_COLORBAR)
figure(gcf);
%im=imagesc(grid.lon, flipud(grid.lat), flipud(reshape(power_vec, grid.Nlat, grid.Nlon)));
imagesc(grid.lon, (grid.lat), (reshape(power_vec, length(grid.lat), length(grid.lon))));
if (nargin<4 || isempty(SHOW_COLORBAR))
    SHOW_COLORBAR = true;
end
if SHOW_COLORBAR
    if (nargin>=3 && ~isempty(tick_vals))
        c=colorbar('Ticks',tick_vals);
    else
        c=colorbar;
    end
    c.Label.String='Power [dBm]';
    c.Label.Interpreter='latex'; 
    c.Label.FontSize = 14;
    c.TickLabelInterpreter = 'latex';
end
xlabel('Lon. [deg]'); ylabel('Lat. [deg]');
set(get(gcf,'Children'),'YDir','normal');
end