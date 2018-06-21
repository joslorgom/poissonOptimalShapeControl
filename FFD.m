function [Nout, P] = FFD(Ni, Nsi, Nsf, Xmin, Xmax, M, N, option)

% Displacement of nodes Ns in plane ST
Ssi = X2S(Nsi, Xmin, Xmax);
Ssf = X2S(Nsf, Xmin, Xmax);
DSs = Ssf - Ssi;

% Initial and final configurations in the plane ST
% figure(100)
% plot(Ssi(:, 1), Ssi(:, 2), 'ro')
% hold on
% plot(Ssf(:, 1), Ssf(:, 2), 'bo')

M = M - 1;
N = N - 1;

% Control points coordinates
[Px, Py] = meshgrid((0:M)/M, (0:N)/N);

% Transformation matrix
F = transMatrix(Ssi, M, N, 'incomplete');

if strcmp(option, 'noExtraRows')
    W = eye(size(F, 1));
    DPx = (F'*(W'*W)*F)\(F'*(W'*W)*DSs(:, 1));
    DPy = (F'*(W'*W)*F)\(F'*(W'*W)*DSs(:, 2));
elseif strcmp(option, 'extraRows')
    n = size(F, 1);
    F = [F; eye((M-1)*(N-1))];
    W = eye(n + (M-1)*(N-1));
    W(1:n, 1:n) = 100*eye(n);
    DPx = (F'*(W'*W)*F)\(F'*(W'*W)*[DSs(:, 1); zeros((M-1)*(N-1), 1)]);
    DPy = (F'*(W'*W)*F)\(F'*(W'*W)*[DSs(:, 2); zeros((M-1)*(N-1), 1)]);
end

k = 1;
for j = 2:M
    for i = 2:N
        Px(i, j) = Px(i, j) + DPx(k);
        Py(i, j) = Py(i, j) + DPy(k);
        k = k + 1;
    end
end

Q = mat2vec(Px, Py);
P = S2X(Q, Xmin, Xmax);

S = X2S(Ni, Xmin, Xmax);
F = transMatrix(S, M, N, 'incomplete');

DS = [F*DPx, F*DPy];

Nout = S2X(S + DS, Xmin, Xmax);

end