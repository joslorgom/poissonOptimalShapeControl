function [Nout] = RBF(nodes, nodes2, I, h, RBFtype, n)

% RBF interpolation
% s = sum_{i=1}^N gamma_i phi( || x - xs_i || ) + p(x)
% xs - coordinates of the N source nodes where displacements are known

% Space dimension
D = 2;
% Degree of polynomial p(x)
deg = 1;

% Coordinates of nodes with prescribed displacements
xs = nodes(I, :);
% Displacements
gs = nodes2(I, :) - xs;
% Number of nodes with prescribed displacements
N = size(xs, 1);

% Matrices initialitation
ncol = 1 + deg*D;
r = zeros(N);
%M = zeros(N);
P = ones(N, ncol);

% Build matrices to find the RBF interpolation parameters
for i = 1:N
    r(i, i) = 0;
    %M(i, i) = radialFunction(0, h, RBFtype, n);
    for j = (i+1):N
        %r = norm( xs(i, :) - xs(j, :) );
        r(i, j) = norm( xs(i, :) - xs(j, :) );
        r(j, i) = r(i, j);
        %M(i, j) = radialFunction(r, h, RBFtype, n);
        %M(j, i) = M(i, j);
    end
    P(i, :) = [1, xs(i, :)];
end
M = radialFunction(r, h, RBFtype, n);

A = [M, P; P', zeros(ncol)];

% Build right-hand side of the system Ax = b
% bx = [gs(:, 1); zeros(ncol, 1)];
% by = [gs(:, 2); zeros(ncol, 1)];
b = [gs; zeros(ncol, 2)];

% Solve for the interpolation parameters
% fParx = A\bx;
% fPary = A\by;
fPar = A\b;

r = zeros(N);
% Build matrices to move all the mesh nodes
for i = 1:size(nodes, 1)
    for j = 1:N
        %r = norm( nodes(i, :) - xs(j, :) );
        r(i, j) = norm( nodes(i, :) - xs(j, :) );
        %M(i, j) = radialFunction(r, h, RBFtype, n);
    end
    P(i, :) = [1, nodes(i, :)];
end
M = radialFunction(r, h, RBFtype, n);

clear A;
A = [M, P];

% Update nodes coordinates
% Nout = nodes + [A*fParx, A*fPary];
Nout = nodes + A*fPar;
        
end