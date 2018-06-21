function [No, I] = inRectangle(Ni, Xmin, Xmax)

% inRectangle Nodes inside the rectangle defined by Xmin and Xmax
% [No, I] = inRectangle(Ni, Xmin, Xmax) returns the nodes and indices in 
% the matrix Ni whose coordinates are inside the rectangle defined by the 
% vectors Xmin and Xmax

a = Ni <= Xmax;
b = Ni >= Xmin;
mask = prod(a.*b, 2);
I = find(mask);
No = Ni(I, :);

end