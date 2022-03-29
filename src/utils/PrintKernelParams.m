function PrintKernelParams(sDistParams,sKernelParams)
dim = sDistParams.dim;
fprintf('*********************************************************\n');
fprintf('* CalcKernelParams\n');
fprintf('*********************************************************\n');
fprintf('omega    (kernel width)            = %8.3f\n', sKernelParams.omega);
fprintf('---------------------------------------------------------\n');
for c = 1:sDistParams.estNumComponents
    fprintf('mu(%2d)    (pdf mean)                 = [',c);
    for d = 1:dim
        fprintf('%8.3f, ', sDistParams.mu{c}(d));
        if (mod(d,8) == 0) && d < dim
            %                 fprintf('\n\t');
            fprintf('...\t'); break;
        end
    end
    fprintf(']\n');
    fprintf('sigma(%2d) (pdf width)                = [',c);
    for d = 1:dim
        fprintf('%8.3f, ', sDistParams.sigma{c}(d));
        if (mod(d,8) == 0) && d < dim
            %                 fprintf('\n\t');
            fprintf('...\t'); break;
        end
    end
    fprintf(']\n');
    fprintf('--> beta(%2d) = 2*sigma(%2d)^2/omega^2 = [', c, c);
    for d = 1:dim
        fprintf('%8.3f, ', sKernelParams.beta{c}(d));
        if (mod(d,8) == 0) && d < dim
            %                 fprintf('\n\t');
            fprintf('...\t'); break;
        end
    end
    fprintf(']\n');
    fprintf('---------------------------------------------------------\n');
end
end