function [Is, K] = springAnalogy(tr, Id)

% Mesh nodes
nodes = tr.Points;
% Number of nodes
N = size(nodes, 1);

% Indices of interior nodes
aux = ones(N, 1);
aux(Id) = 0;
Is = find(aux);

% Elements
elements = tr.ConnectivityList;
% Number of elements
Ne = size(elements, 1);

% Stiffness matrix
K = sparse(N, N);

% Run over all mesh elements
for i = 1:Ne
    % Nodes in element i
    n = zeros(size(elements, 2), size(nodes, 2));
    for j = 1:size(elements, 2)
        n(j, :) = nodes(elements(i, j), :);
    end
    l12 = norm(n(1, :) - n(2, :));
    l13 = norm(n(1, :) - n(3, :));
    l23 = norm(n(2, :) - n(3, :));
    
    n1 = elements(i, 1);
    n2 = elements(i, 2);
    n3 = elements(i, 3);
    
    K(n1, n2) = -1/l12;
    K(n2, n1) = -1/l12;
    K(n1, n3) = -1/l13;
    K(n3, n1) = -1/l13;
    K(n2, n3) = -1/l23;
    K(n3, n2) = -1/l23;
    
%     K(n1, n1) = K(n1, n1) + 1/l12 + 1/l13;
%     K(n2, n2) = K(n2, n2) + 1/l12 + 1/l23;
%     K(n3, n3) = K(n3, n3) + 1/l13 + 1/l23;
end

for i = 1:N
    K(i, i) = -sum(K(i, :));
end

end