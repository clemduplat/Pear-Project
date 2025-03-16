function value = partFvAddition(triangleAir, r1, r2, r3)
    matrix = [r3 - r1; r1 + 2*r2 + r3; r1 + r2 + 2*r3];
    value = triangleAir*matrix/12;
end