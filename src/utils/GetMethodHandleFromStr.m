function funcHandle = GetMethodHandleFromStr(methodStr)
switch methodStr
    case 'RoMix'
        funcName = 'InterpGraphSignalRoMix';
    case 'Rep. Thm.'
        funcName = 'InterpGraphSignalRepThm';
    case 'VSPW'
        funcName = 'InterpGraphSignalPesenson';
    case 'Nystr\"{o}m'
        funcName = 'InterpGraphSignalNystrom';
    case 'w-kNN'
        funcName = 'InterpGraphSignalKnn';
    otherwise
        error('invalid methodStr')
end
funcHandle = str2func(funcName);
end