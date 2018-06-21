function [dist] = distance(nodes, Ic)

% distance returns the minimum distance of all nodes to the boundary
% described by the nodes in Ic

M = length(Ic);
N = size(nodes, 1);

aux = zeros(N, M);

for j = 1:M
    nc = nodes(Ic(j), :);
    dc = nodes - nc;
    aux(:, j) = sqrt(dc(:, 1).^2 + dc(:, 2).^2);
end

dist = min(aux, [], 2);

end