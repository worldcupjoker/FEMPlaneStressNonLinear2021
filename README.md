FEMPlaneStressNonLinear2021

Main file: NonLinearQ8.m

This file solves a non linear 2-D plane stress problem.

Create an object: obj = NonLinearQ8('Prob2_133_big.inp', nu, E0, q, mode)

u_xx = obj.getUxx();
u_yy = obj.getUyy();

A main function file is provided to perform sample calculation.
