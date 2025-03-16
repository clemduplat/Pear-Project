function value = H2(length, r1, r2)
    matrix = [2*r1 + r2; r1 + 2*r2];
    value = length*matrix/6;
end