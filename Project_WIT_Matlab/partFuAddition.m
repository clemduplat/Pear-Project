function value = partFuAddition(triangleAir, r1, r2, r3)
    matrix = [[(6*r1 + 2*r2 + 2*r3) (2*r1 + 2*r2 + r3) (2*r1 + r2 + 2*r3)]; [(2*r1 + 2*r2 + r3) (2*r1 + 6*r2 + 2*r3) (r1 + 2*r2 + 2*r3)] ;[(2*r1 + r2 + 2*r3) (r1 + 2*r2 + 2*r3) (2*r1 + 2*r2 + 6*r3)]];
    value = triangleAir*matrix/60;
end
