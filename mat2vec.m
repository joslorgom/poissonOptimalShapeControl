function [P] = mat2vec(Px, Py)

[N, M] = size(Px);

P = zeros(M*N, 2);
k = 1;
for j = 1:M
    for i = 1:N
        P(k, :) = [Px(i, j), Py(i, j)];
        k = k + 1;
    end
end

end