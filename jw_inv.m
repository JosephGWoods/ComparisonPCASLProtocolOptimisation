function inverse = jw_inv(m)

% Analytical inversion of 2x2 matrix

% Calculate the determinant
d = (m(1)*m(4)) - (m(2)*m(3));

% Calculate the inverse
inverse = (1/d) .* [m(4), -m(2); -m(3), m(1)];

end