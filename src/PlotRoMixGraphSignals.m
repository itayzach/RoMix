function PlotRoMixGraphSignals(sPlotParams, sPreset, sDataset, sKernelParams, mSigCnvrtRec, mSigCnvrtInt, C)
if sPlotParams.b_globalPlotEnable
    mSigCnvrtRecRef = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);
    mSigCnvrtIntRef = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.yt);
    if sPreset.b_runGraphSignals
        PlotGraphSignalsWrapper(sPlotParams, sPreset, sKernelParams, sDataset, ...
            mSigCnvrtRecRef, mSigCnvrtIntRef, mSigCnvrtRec, mSigCnvrtInt, C, 'RoMix')
    end
    if sPreset.b_interpEigenvecs
        PlotEigenvectorsWrapper(sPlotParams, sPreset, sDataset, mSigCnvrtRecRef, mSigCnvrtIntRef, mSigCnvrtRec, mSigCnvrtInt);
    end
end
end