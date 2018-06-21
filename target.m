function [ud] = target(X)

x = X(1);
y = X(2);

u1 = 0;
u2 = 1;
R1 = 0.3;
R2 = 1;
f = 0;

r = sqrt(x^2 + y^2);

if r <= R1
    ud = u1;
else
    c1 = ( u2 - u1 + f/4*( R2^2 - R1^2 ) ) / log(R2/R1);
    c2 = ( u1 + u2 )/2 + f/8*( R1^2 + R2^2 ) - c1/2*log(R1*R2);
    ud = -f/4*r^2 + c1*log(r) + c2;
end

% if x >= 3/4 || x <= 1/4 || y >= 3/4 || y <= 1/4
%     ud = sin(4*pi*x)*sin(4*pi*y);
% else
%     ud = 0;
% end

end