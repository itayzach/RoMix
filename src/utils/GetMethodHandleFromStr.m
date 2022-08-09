function funcHandle = GetMethodHandleFromStr(methodStr)
switch methodStr
    case 'RoMix'
        funcName = 'InterpGraphSignalRoMix';
    case 'Rep. Thm.'
        funcName = 'InterpGraphSignalRepThm';
    case 'PW'
        funcName = 'InterpGraphSignalPesenson';
    case 'Nystrom'
        funcName = 'InterpGraphSignalNystrom';
    case 'kNN'
        funcName = 'InterpGraphSignalKnn';
    otherwise
        error('invalid methodStr')
end
funcHandle = str2func(funcName);
end