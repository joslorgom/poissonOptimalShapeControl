function [X] = S2X(S, Xmin, Xmax)

% Transform the nodes coordinates in plane ST into the plane XY
X = Xmin + (Xmax - Xmin).*S;

end