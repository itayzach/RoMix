function sResults = SaveResults(sPreset, sPlotParams, mAccRec, mAccStdRec, mAccInt, mAccStdInt, vTrainTime, vTrainTimeStd, vIntTime, vIntTimeStd)
sResults.sPreset         = sPreset;
sResults.sPlotParams     = sPlotParams;
%sResults.tSigCnvrtRec    = tSigCnvrtRec;
%sResults.tSigCnvrtInt    = tSigCnvrtInt;
%sResults.tSigCnvrtRecRef = tSigCnvrtRecRef;
%sResults.tSigCnvrtIntRef = tSigCnvrtIntRef;
sResults.mAccRec         = mAccRec;
sResults.mAccStdRec      = mAccStdRec;
sResults.mAccStdInt      = mAccStdInt;
sResults.mAccInt         = mAccInt;
sResults.vTrainTime      = vTrainTime;
sResults.vTrainTimeStd   = vTrainTimeStd;
sResults.vIntTime        = vIntTime;
sResults.vIntTimeStd     = vIntTimeStd;
end