clc
clear

A = xlsread('\test.xlsx');
x = A(:, 1);
y = A(:, 2);
z = A(:, 3);
scatter3(x,y,z)