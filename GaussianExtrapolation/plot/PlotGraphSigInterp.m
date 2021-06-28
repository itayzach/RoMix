function PlotGraphSigInterp(sSimParams, dataDim, G, samplingRatio, sampleInd, f, gsp_interp_signal, f_int_no_transform, graphName)
fig = figure;
param.show_edges = false;
subplot(131);
    if dataDim > 1
        gsp_plot_signal(G,f,param); 
    else
        plot(G.coords, f, '.');
    end
    if sSimParams.b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G.coords(sampleInd), f(sampleInd), 'ko')
        elseif dataDim == 2
            scatter(G.coords(sampleInd,1), G.coords(sampleInd,2), 'ko')
        elseif dataDim == 3
            scatter3(G.coords(sampleInd,1), G.coords(sampleInd,2), G.coords(sampleInd,3), 'ko')
        end
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Original graph-signal' newline ...
           '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f') newline ...
           'Sampling ratio = ' num2str(samplingRatio, '%.2f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(132); 
    if dataDim > 1
        gsp_plot_signal(G,gsp_interp_signal,param); 
    else
        plot(G.coords, gsp_interp_signal, '.');
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Pesenson''s interpolation' newline ...
        '$f_{{\bf int}}^T L f_{{\bf int}}$ = ' num2str(gsp_interp_signal'*G.L*gsp_interp_signal, '%.3f') newline ...
        'error: ' num2str(norm(f - gsp_interp_signal)/norm(f), '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(133); 
    if dataDim > 1
        gsp_plot_signal(G,f_int_no_transform,param);
    else
        plot(G.coords, f_int_no_transform, '.');
    end
    set(gca,'FontSize', 14);
    title(['Our interpolation (no transformation)' newline '$f_{{\bf int}}^T L f_{{\bf int}}$ = ' num2str(f_int_no_transform'*G.L*f_int_no_transform, '%.3f') newline ...
        'error: ' num2str(norm(f - f_int_no_transform)/norm(f), '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
    view(0,90)
set(gcf,'Position', [400 100 1200 400])   
if isfield(sSimParams, 'outputFolder')
    saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_interp_no_transform_', num2str(round(samplingRatio*100), '%d'), 'prec_sampling'), 'epsc');
end

end

