function [M, K, Cx, Cy, Ae, An, fse] = assembleMatrices(tr)

elem = tr.ConnectivityList;
nodes = tr.Points;

Ne = size(elem, 1);
Nn = size(nodes, 1);

% Stiffness Matrix
K = sparse(Nn, Nn);
% Mass Matrix
M = sparse(Nn, Nn);
% Convection Matrices
Cx = sparse(Nn, Nn);
Cy = sparse(Nn, Nn);

Ae = zeros(Ne, 1);
An = zeros(Nn, 1);

fse = zeros(Ne, 1);

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
    
    Cxe = zeros(3);
    Cxe(:, 1) = -(b(2) + b(3))/6;
    Cxe(:, 2) = b(2)/6;
    Cxe(:, 3) = b(3)/6;
    
    Cye = zeros(3);
    Cye(:, 1) = -(c(2) + c(3))/6;
    Cye(:, 2) = c(2)/6;
    Cye(:, 3) = c(3)/6;
    
    for i = 1:3
        %F(n(i), 1) = f(x(i, :));
        for j = 1:3
            K(n(i), n(j)) = K(n(i), n(j)) + (b(i)*b(j) + c(i)*c(j))/(4*Ae(k));
            M(n(i), n(j)) = M(n(i), n(j)) + ((i == j) + 1)/24*Ae(k)*2;
            Cx(n(i), n(j)) = Cx(n(i), n(j)) + Cxe(i, j);
            Cy(n(i), n(j)) = Cy(n(i), n(j)) + Cye(i, j);
        end
    end
    
    An(n) = An(n) + 1/Ae(k);
    
    Ak = [(x(2, :) - x(1, :))', (x(3, :) - x(1, :))'];
    Lk = Ak'*Ak;
    alpha = sqrt( Lk(1, 1)*Lk(2, 2) - Lk(1, 2)^2 );
    fse(k) = sqrt(3)*alpha/( Lk(1, 1) + Lk(2, 2) - Lk(1, 2) );
    
end

end