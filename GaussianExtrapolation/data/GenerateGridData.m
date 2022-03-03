function data = GenerateGridData(dim, Nx, xMin, xMax)
if dim == 1
    dx = (xMax-xMin)/Nx;
    data = ((1:Nx)*dx - dx/2).';
elseif dim == 2
    dx = (xMax-xMin)./Nx;
    x = ((1:Nx(1))*dx(1) - dx(1)/2).';
    y = ((1:Nx(2))*dx(2) - dx(2)/2).';
    [X, Y] = meshgrid(x, y);
    data = [X(:) Y(:)];
elseif dim == 3
    dx = (xMax-xMin)./Nx;
    x = ((1:Nx(1))*dx(1) - dx(1)/2).';
    y = ((1:Nx(2))*dx(2) - dx(2)/2).';
    z = ((1:Nx(3))*dx(3) - dx(3)/2).';
    [X, Y, Z] = meshgrid(x, y, z);
    data = [X(:) Y(:) Z(:)];
else
    error('not supported');
end
end