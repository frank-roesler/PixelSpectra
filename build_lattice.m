function L = build_lattice(z1, z2, h)
%   Defines lattice L of distance h in the complex plane given the corner 
%   points z1, z2 as input.
    r1 = real(z1);
    r2 = real(z2);
    i1 = imag(z1);
    i2 = imag(z2);
    X = r1:h:r2;
    Y = i1:h:i2;
    [XX,YY] = meshgrid(X,Y);
    L = XX+1i*YY;
end