function PlotNumEigsPerComp(sKernelParams)
nEigs = numel(sKernelParams.vComponentIndex);
figure; histogram(sKernelParams.vComponentIndex);
hold on; plot(sKernelParams.sDistParams.GMModel.ComponentProportion*nEigs)
legend('Hist', 'CompProp')
title('\# of eigs in each component')
end