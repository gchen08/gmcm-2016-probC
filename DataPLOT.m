clc
clear

A = xlsread('C:\Users\Gong\workspace\ProblemC\data.xlsx');
y = A(:, 1);
x = 1:1:length(y);
plot(x, y, '+b');
hold on
z = mean(y)
plot(x, z, 'r','LineWidth',8);

