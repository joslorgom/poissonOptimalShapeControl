function [F] = assembleVector(f, tr)

nodes = tr.Points;
elem = tr.ConnectivityList;

Nn = size(nodes, 1);
Ne = size(elem, 1);

% F = zeros(Nn, 1);
% for i = 1:Nn
%     F(i, 1) = f(nodes(i, :));
% end

F = zeros(Nn, 1);

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
    
    xq = [x(1, :)/2 + x(3, :)/2;
        x(1, :)/2 + x(2, :)/2;
        x(2, :)/2 + x(3, :)/2];
    
    for i = 1:3
        F(n(i), 1) = F(n(i), 1) + f(xq(i, :))/6*Ae(k)*2;
    end
end

end