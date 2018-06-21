function [K, F] = dirichletBC(K, F, Ic, Ie, uc, ue)

K(Ic, :) = 0;
K(Ic, Ic) = eye(size(Ic, 1));
F(Ic) = uc;

K(Ie, :) = 0;
K(Ie, Ie) = eye(size(Ie, 1));
F(Ie) = ue;

end