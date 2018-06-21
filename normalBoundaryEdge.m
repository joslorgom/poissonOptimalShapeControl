function [nbe, xbe] = normalBoundaryEdge(tr, fb)

nodes = tr.Points;

Nbe = size(fb, 1);

xbe = zeros(Nbe, 2);
nbe = zeros(Nbe, 2);

for i = 1:Nbe
    i1 = fb(i, 1);
    i2 = fb(i, 2);
    xbe(i, :) = 0.5 * ( nodes(i1, :) + nodes(i2, :) );
    dn = nodes(i2, :) - nodes(i1, :);
    nbe(i, :) = [dn(2), -dn(1)]/norm(dn);
end

% figure(2)
% hold on
% plot(xbe(:, 1), xbe(:, 2), 'ro')
% quiver(xbe(:, 1), xbe(:, 2), nbe(:, 1), nbe(:, 2), 'r', 'AutoScaleFactor', 0.2)

end