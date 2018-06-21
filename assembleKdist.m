function [K] = assembleKdist(tr, gamma)

elem = tr.ConnectivityList;
nodes = tr.Points;

Ne = size(elem, 1);
Nn = size(nodes, 1);

% Stiffness Matrix
K = sparse(Nn, Nn);

Ae = zeros(Ne, 1);

for k = 1:Ne
    
    n = [elem(k, 1);
        elem(k, 2);
        elem(k, 3)];
    
    x = [nodes(n(1), :);
        nodes(n(2), :);
        nodes(n(3), :)];
    
    a = [x(2, :) - x(3, :);
        x(3, :) - x(1, :);
        x(1, :) - x(2, :)];
    
    Ae(k) = (-a(3, 1)*a(2, 2) + a(2, 1)*a(3, 2))/2;
    
    b = a(:, 2);
    c = -a(:, 1);
    
    D = ( gamma(n(1)) + gamma(n(2)) + gamma(n(3)) )/3;
    
    for i = 1:3
        for j = 1:3
            K(n(i), n(j)) = K(n(i), n(j)) + D * (b(i)*b(j) + c(i)*c(j))/(4*Ae(k));
        end
    end
    
end

end