function value = KvAddition(triangleAir, sig_vr, sig_vz, r1, r2, r3, z1, z2, z3 )
    matrix1 = [[(z2 - z3) (r3 - r2)]; [(z3 - z1) (r1 - r3)] ;[(z1 - z2) (r2 - r1)]];
    matrix2 = diag([sig_vr, sig_vz]);
    value = (r1 + r2 + r3)*matrix1*matrix2*(matrix1')/(12*triangleAir);
end