function [nbv, xbv] = normalBoundaryVertex(tr, fb, nbe)

nodes = tr.Points;

Nbe = size(fb, 1);

xbv = zeros(Nbe, 2);
nbv = zeros(Nbe, 2);
%fv = zeros(Nfb, 2);
for i = 1:Nbe
    k1 = find(fb(:, 1) == i);
    %i1 = fb(k, 2);
    k2 = find(fb(:, 2) == i);
    %i2 = fb(k, 1);
    %fv(i, :) = [k1, k2];
    nsum = nbe(k1, :) + nbe(k2, :);
    nbv(i, :) = nsum/norm(nsum);
    xbv(i, :) = nodes(i, :);
end

% figure(2)
% hold on
% plot(xbv(:, 1), xbv(:, 2), 'bo')
% quiver(xbv(:, 1), xbv(:, 2), nbv(:, 1), nbv(:, 2), 'b', 'AutoScaleFactor', 0.2)
    

end