function data = GenerateGridData(dim, n, xMin, xMax)
if dim == 1
    data = linspace(xMin(1),xMax(1),n)';
elseif dim == 2
    x = linspace(xMin(1),xMax(1),round(sqrt(n)));
    y = linspace(xMin(2),xMax(2),round(sqrt(n)));
    [X, Y] = meshgrid(x, y);
    data = [X(:) Y(:)];
else
    error('not supported');
end
end