function [] = plot_classifier(classifier, output_folder, xTest, yTest)

x0     = 10;
y0     = 250;
width  = 550;
height = 700;
if strcmp(classifier.Name, 'laprlsc') || strcmp(classifier.Name, 'clustering')
    %--------------------------------------------------------------------------
    % Plot LapRLS
    %--------------------------------------------------------------------------
    alpha_laprls = classifier.alpha;
    xTrain_laprls = classifier.xtrain;
    yTrain_laprls = classifier.ytrain;
    xMax = max(max(xTrain_laprls,[],1));
    xMin = min(min(xTrain_laprls,[],1));
    step = (xMax - xMin)/100;
    x1 = xMin:step:xMax;
    x2 = x1;
    [XX1,XX2] = meshgrid(x1,x2);

    X=[XX1(:) XX2(:)];
    tic;
    mK_xTrain_X = calckernel(classifier.Kernel,classifier.KernelParam, xTrain_laprls, X);
    run_time = toc;
    fprintf('K time = %f\n', run_time);
    mK_xTrain_xTrain = calckernel(classifier.Kernel,classifier.KernelParam, xTrain_laprls, xTrain_laprls);
    vKa = mK_xTrain_X*alpha_laprls;
    vKa_train = mK_xTrain_xTrain*alpha_laprls;
    mKa = reshape(vKa,length(x1),length(x2));
    
    fig = figure;
    sgtitle(['LapRLS with $n$ = ' num2str(length(classifier.alpha)) newline ...
        '$\gamma_A = ' num2str(classifier.gammas(1), '%.4f') '$ $\gamma_I = ' num2str(classifier.gammas(2), '%.4f') '$' newline ...
        'Test error = ' num2str(classifier.test_error, '%.1f'), '$\%$' ], ...
        'Interpreter', 'latex', 'FontSize', 14);
    subplot(2,1,1)
    surf(XX1,XX2,mKa, 'edgecolor', 'none')
    hold on;
    pos=find(yTrain_laprls==1);
    neg=find(yTrain_laprls==-1);
    unlab=find(yTrain_laprls==0);
    scatter3(xTrain_laprls(unlab,1), xTrain_laprls(unlab,2), vKa_train(unlab), 'ks');
    scatter3(xTrain_laprls(pos,1),xTrain_laprls(pos,2), vKa_train(pos),200, 'rd','MarkerFaceColor','r','MarkerEdgeColor','k'); hold on;
    scatter3(xTrain_laprls(neg,1),xTrain_laprls(neg,2), vKa_train(neg),200, 'bo' ,'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;
    
    xlabel('$x_1$', 'Interpreter', 'latex')
    ylabel('$x_2$', 'Interpreter', 'latex')
    zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
    colorbar;
    set(gca,'FontSize', 14);
    subplot(2,1,2)
    contourf(XX1,XX2,mKa,[0 0]);shading flat;
    hold on;
    plot2D(xTrain_laprls,yTrain_laprls,15,'ks');
    plot2D(xTest,yTest,15,'k*');
    set(gca,'FontSize', 14);
    set(gcf,'Position',[x0+width y0 width height])
    saveas(fig,[output_folder filesep 'fig_' classifier.Name], 'epsc');
    
elseif strcmp(classifier.Name, 'eigrls')
    %--------------------------------------------------------------------------
    % Plot EigRLS (Phi*c)
    %--------------------------------------------------------------------------
    c = classifier.c;
    sParams = classifier.sParams;
    xTrain_eigrls = classifier.xtrain;
    yTrain_eigrls = classifier.ytrain;
    xMax = max(max(xTrain_eigrls,[],1));
    xMin = min(min(xTrain_eigrls,[],1));
    step = (xMax - xMin)/100;
    x1 = xMin:step:xMax;
    x2 = x1;
    [XX1,XX2] = meshgrid(x1,x2);

    X=[XX1(:) XX2(:)];    
    mPhi_m_X = zeros(length(X), sParams.ExtrplM-sParams.FirstM);
    
    tic;
    for i = sParams.FirstM:sParams.ExtrplM-1
        m = OneDim2TwoDimIndex(sParams.multindexToSingleIndexMap(i+1)-1);
        mPhi_m_X(:, i-sParams.FirstM+1) = phi(sParams,m,X);
    end
    run_time = toc;
    fprintf('Phi time = %f\n', run_time);
    
%     mPhi_m_xTrain = zeros(length(xTrain_eigrls), sParams.ExtrplM-sParams.FirstM);    
% %     mPhi_m_xTest = zeros(length(xTest),sParams.ExtrplM-sParams.FirstM);
% %     vLambda = zeros(1, sParams.ExtrplM-sParams.FirstM);
%     for i = sParams.FirstM:sParams.ExtrplM-1
%         m = OneDim2TwoDimIndex(sParams.multindexToSingleIndexMap(i+1)-1);
% %         vLambda(i-sParams.FirstM+1) = lambda(sParams,m);
%         mPhi_m_xTrain(:, i-sParams.FirstM+1) = phi(sParams,m,xTrain_eigrls);
% %         mPhi_m_xTest(:, i+1) = phi(sParams,m,xTest);
%     end
    
    mPhi_m_xTrain = classifier.mPhi_m_xTrain;
    
    vPhi_X_c = mPhi_m_X*c;
    vPhi_xTrain_c = mPhi_m_xTrain*c;
    mPhi_X_c = reshape(vPhi_X_c,length(x1),length(x2));
    % mPhi_xt_c = reshape(vPhi_xt_c,length(x1),length(x2));
    
    fig = figure;
    sgtitle(['EigRLS with $M$ = ' num2str(sParams.ExtrplM-sParams.FirstM) newline ...
        '$\omega = $ ' num2str(sParams.omega) ';   $\gamma_A = ' num2str(classifier.gammas(1), '%.4f') '$ $\gamma_I = ' num2str(classifier.gammas(2), '%.4f') '$' newline ...
        'Test error = ' num2str(classifier.test_error, '%.1f'), '$\%$' ], ...
        'Interpreter', 'latex', 'FontSize', 14);
    subplot(2,1,1)
    surf(XX1,XX2,mPhi_X_c, 'edgecolor', 'none')
    hold on;
    pos=find(yTrain_eigrls==1);
    neg=find(yTrain_eigrls==-1);
    unlab=find(yTrain_eigrls==0);
    scatter3(xTrain_eigrls(unlab,1), xTrain_eigrls(unlab,2), vPhi_xTrain_c(unlab), 'ks');
    scatter3(xTrain_eigrls(pos,1),xTrain_eigrls(pos,2), vPhi_xTrain_c(pos),200, 'rd','MarkerFaceColor','r','MarkerEdgeColor','k'); hold on;
    scatter3(xTrain_eigrls(neg,1),xTrain_eigrls(neg,2), vPhi_xTrain_c(neg),200, 'bo' ,'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;
    xlabel('$x_1$', 'Interpreter', 'latex')
    ylabel('$x_2$', 'Interpreter', 'latex')
    zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
    colorbar;
    set(gca,'FontSize', 14);
    subplot(2,1,2)
    contourf(XX1,XX2,mPhi_X_c,[0 0]);shading flat;
    hold on;
    plot2D(xTrain_eigrls,yTrain_eigrls,15,'ks');
    plot2D(xTest,yTest,15,'k*');
    set(gca,'FontSize', 14);
    set(gcf,'Position',[x0+width+width y0 width height])
    saveas(fig,[output_folder filesep 'fig_eigrls_M_' num2str(sParams.ExtrplM)], 'epsc');
else
    error('unknown classifier name')
end
end