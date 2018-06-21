function [Cc2, wtan] = tangential(nodes, Cc)

% Tangential receives the mesh nodes and a sorted list with the indices of
% the mesh nodes lying on the controlled boundary, whose first and second 
% indices coincide. It returns a shortened list and the displacements to be
% applied to these nodes in order to preserve the shape of the elements
% having any edge lying on the controlled boundary.

wrapN = @(x, n) (1 + mod(x-1, n));

% Shorten list
Cc2 = Cc(1:end-1);
M = size(Cc2, 1);

% Edges lengths
Lc = zeros(M, 1);
% Tangent vectors
tc = zeros(M, size(nodes, 2));

for j = 1:M
    Lc(j) = norm(nodes(Cc(j+1), :) - nodes(Cc(j), :));
    tc(j, :) = nodes(Cc2(wrapN(j+1, M)), :) - nodes(Cc2(wrapN(j-1, M)), :);
    tc(j, :) = tc(j, :)/norm(tc(j, :));
end
% Average edge length
Lav = sum(Lc)/M;

% Compute corrector displacement
dtau = zeros(M, 1);
wtan = zeros(M, 2);
for j = 2:M
    dtau(j) = dtau(j-1) + (Lav - Lc(j-1));
    wtan(j, :) = dtau(j) * tc(j, :);
end

end