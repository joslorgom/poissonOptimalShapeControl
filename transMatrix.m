function [F] = transMatrix(S, M, N, option)

% Number of nodes to be displaced
L = size(S, 1);

% Evaluate Bernstein polynomials in S = [s, t]
% bi = (M i) s^i (1-s)^(M-i); i = 0,...,M
Bi = zeros(L, M+1);
for i = 0:M
    Bi(:, i+1) = nchoosek(M, i) * S(:, 1).^i .* (1 - S(:, 1)).^(M - i);
end
% bj = (N j) t^j (1-t)^(N-j); j = 0,...,N
Bj = zeros(L, N+1);
for j = 0:N
    Bj(:, j+1) = nchoosek(N, j) * S(:, 2).^j .* (1 - S(:, 2)).^(N - j);
end

if strcmp(option, 'complete')

    F = zeros(L, (M+1)*(N+1));
    k = 1;
    for i = 1:(M+1)
        for j = 1:(N+1)
            F(:, k) = Bi(:, i) .* Bj(:, j);
            k = k + 1;
        end
    end
    
elseif strcmp(option, 'incomplete')

    F = zeros(L, (M-1)*(N-1));
    k = 1;
    for i = 2:M
        for j = 2:N
            F(:, k) = Bi(:, i) .* Bj(:, j);
            k = k + 1;
        end
    end
    
else
    F = transMatrix(S, M, N, 'complete');    
end

end