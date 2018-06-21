function [Cc] = sortNodes(fb, Ic)

Nc = size(Ic, 1);
Cc = zeros(Nc+1, 1);
Cc(1) = Ic(1);
v = Cc(1);
for i = 1:Nc
    k = find(fb(:, 1) == v);
    v = fb(k, 2);
    Cc(i + 1) = v;
end

end