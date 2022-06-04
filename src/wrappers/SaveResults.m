function sResults = SaveResults(sPreset, sPlotParams, tSigCnvrtRec, tSigCnvrtInt, tSigCnvrtRecRef, tSigCnvrtIntRef, mAccRec, mAccStdRec, mAccInt, mAccStdInt, tTrainTime, tIntTime)
sResults.sPreset         = sPreset;
sResults.sPlotParams     = sPlotParams;
sResults.tSigCnvrtRec    = tSigCnvrtRec;
sResults.tSigCnvrtInt    = tSigCnvrtInt;
sResults.tSigCnvrtRecRef = tSigCnvrtRecRef;
sResults.tSigCnvrtIntRef = tSigCnvrtIntRef;
sResults.mAccRec         = mAccRec;
sResults.mAccStdRec      = mAccStdRec;
sResults.mAccStdInt      = mAccStdInt;
sResults.mAccInt         = mAccInt;
sResults.tTrainTime      = tTrainTime;
sResults.tIntTime        = tIntTime;
end