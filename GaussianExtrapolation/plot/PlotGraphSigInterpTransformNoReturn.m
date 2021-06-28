function PlotGraphSigInterpTransformNoReturn(sSimParams, dataDim, G, f, G_tilde, f_tilde_interp, samplingRatio, sample_ind, graphName, verticesTransformation, interpMethod)
fig = figure;
param.show_edges = false;
subplot(1,2,1);
    if dataDim > 1
        gsp_plot_signal(G,f,param); 
    else
        plot(G.coords, f, '.');
    end
    if sSimParams.b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G.coords(sample_ind), f(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Original graph-signal' newline ...
           '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f') newline...
           'Sampling ratio = ' num2str(samplingRatio, '%.2f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(1,2,2); 
    if dataDim > 1
        gsp_plot_signal(G_tilde,f_tilde_interp,param); 
    else
        plot(G_tilde.coords, f_tilde_interp, '.');
        xlabel('$\tilde{v}$', 'Interpreter', 'latex', 'FontSize', 14);
    end
    if sSimParams.b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G_tilde.coords(sample_ind), f_tilde_sampled_padded(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G_tilde.coords(sample_ind,1), G_tilde.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G_tilde.coords(sample_ind,1), G_tilde.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end    
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Interpolated graph-signal on $\tilde{G}$' newline ...
           'using ' interpMethod newline...
           '$\tilde{f}_{\bf int}^T \tilde{L} \tilde{f}_{\bf int}$ = ' num2str(f_tilde_interp'*G_tilde.L*f_tilde_interp, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
set(gcf,'Position', [400 100 800 400])   
if isfield(sSimParams, 'outputFolder')
    saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_', verticesTransformation, '_', interpMethod, '_interp_with_transform_no_return_', num2str(round(samplingRatio*100), '%d'), 'prec_sampling'), 'epsc');
end
end

