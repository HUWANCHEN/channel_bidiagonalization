function [v,b,ephi] = house(x)
% The first entry of x*ephi is real-valued
% ephi*x-b*v*v'*ephi*x = [||x||;0;0;...;0]

if size(x,2) > 1
    error('x should be a column vector !');
end
n = length(x);
ephi = exp(-1j*angle(x(1)));
x = ephi*x;
if n > 1
    sigm = x(2:n)'*x(2:n);
    v = [1;x(2:n)];    
    if sigm == 0
        b = 0;
    else
        m = norm(x);
        v(1) = -sigm/(x(1)+m);
        b = 2*abs(v(1))^2/(sigm+v(1)^2);
        v = v/v(1);

    end
else
    v = 1;
    b = 0;
end