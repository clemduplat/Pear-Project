function value = H1(length, r1, r2)
    matrix = [[(3*r1 + r2) (r1 + r2)]; [(r1 + r2) (r1 + 3*r2)]];
    value = length*matrix/12;
end