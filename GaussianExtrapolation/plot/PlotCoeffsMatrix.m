function PlotCoeffsMatrix(C1, C1Title, C2, C2Title, C3, C3Title)

if exist('C2', 'var')
    cMin = min([C1(:); C2(:)]);
    cMax = max([C1(:); C2(:)]);
else
    cMin = min(C1(:));
    cMax = max(C1(:));
end

figure('Name', 'C');
tiledlayout('flow')
nexttile;
imagesc(C1);
colormap('jet');
colorbar('TickLabelInterpreter', 'latex');
caxis([cMin, cMax]);
title(C1Title, 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

if exist('C2', 'var') && ~isempty(C2)
    nexttile;
    imagesc(C2);
    colormap('jet');
    colorbar('TickLabelInterpreter', 'latex');
    caxis([cMin, cMax]);
    title(C2Title, 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
    
    nexttile;
    imagesc(abs(C1-C2));
    colormap('jet');
    colorbar('TickLabelInterpreter', 'latex');
    title('$|C_1 - C_2|$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
    
    if exist('C3', 'var')
        nexttile;
        imagesc(C3);
        colormap('jet');
        colorbar('TickLabelInterpreter', 'latex');
        title(C3Title, 'interpreter', 'latex', 'FontSize', 16);
        set(gca,'FontSize', 14);
        
    end
end
end