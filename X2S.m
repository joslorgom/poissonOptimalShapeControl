function [S] = X2S(X, Xmin, Xmax)

% Transform the nodes coordinates in plane XY into the square [0, 1]x[0, 1]
S = (X - Xmin)./(Xmax - Xmin);

end