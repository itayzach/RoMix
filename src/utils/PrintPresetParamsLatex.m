function PrintPresetParamsLatex(sPreset)
s = [];
%s = [s, '\begin{table}[htbp]', newline ];
s = [s, '\begin{minipage}[c]{0.5\textwidth}', newline ];
s = [s, '    \centering', newline ];
s = [s, '    \begin{tabular}{ |l|l| } ', newline ];
s = [s, '    \hline', newline ];
s = [s, '    Parameter    & Value  \\ ', newline ];
s = [s, '    \hline\hline', newline ];
s = [s, '    Kernel width & $\omega = ' num2str(sPreset.omegaTilde) '$ \\', newline ];
s = [s, '    \hline', newline ];
s = [s, '    \# GMM comp.    & $k = ' num2str(sPreset.gmmNumComponents) '$  \\ ', newline ];
s = [s, '    \hline', newline ];
if all(floor(sPreset.MTilde./sPreset.gmmNumComponents) == sPreset.MTilde./sPreset.gmmNumComponents)
    s = [s, '    \# eigs. & $M = ' num2str(sPreset.MTilde./sPreset.gmmNumComponents(1)) '\cdot k$  \\ ', newline ];
else
    s = [s, '    \# eigs. & $M = ' num2str(sPreset.MTilde) '$  \\ ', newline ];
end
s = [s, '    \hline', newline ];
s = [s, '    $\| \cdot \|_{\H_K}$ penalty & $\gamma = ' num2str(sPreset.gamma1, '%g') '$\\', newline ];
s = [s, '    \hline', newline ];
if isfield(sPreset.sDatasetParams, 'b_runVAE') && sPreset.sDatasetParams.b_runVAE
    s = [s, '    Latent dim & $d = ' num2str(sPreset.sDatasetParams.latentDim) '$ \\', newline ];
    s = [s, '    \hline', newline ];
end
if numel(sPreset.nLabeled) > 1
    s = [s, '    $|\V_\ell|$ & $\ell = $ varying percentage of $n$\\', newline ];
else
    s = [s, '    $|\V_\ell|$ & $\ell = ' num2str(sPreset.nLabeled, '%d') ' $\\', newline ];
end
s = [s, '    \hline', newline ];
s = [s, '    $|\V|$ & $n = ' num2str(sPreset.n, '%d') ' $\\', newline ];
s = [s, '    \hline', newline ];
s = [s, '    $|\V_N|$ & $N = ' num2str(sPreset.N, '%d') ' $\\', newline ];
s = [s, '    \hline', newline ];

s = [s, '    \end{tabular}', newline ];
%s = [s, '    \captionof{table}{Parameters for ' AddSpaceBeforeUpper(sPreset.verticesPDF) ' dataset}', newline ];
%s = [s, '    \label{tab:params_' sPreset.verticesPDF '}', newline ];
s = [s, '\end{minipage}'];
%s = [s, '\end{table}'];
fprintf('%s\n', s);

end

function spacedStr = AddSpaceBeforeUpper(str)
isUpperCase = isstrprop(str, 'upper');
indVec = 1:numel(str);
if sum(isUpperCase) == numel(str) || sum(isUpperCase) == 1
    spacedStr = str;
else
    upperInd = indVec(isUpperCase);
    spacedStr = [str(1:upperInd(2)-1), ' ', str(upperInd(2):end)];
end

end