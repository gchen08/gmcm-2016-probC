clc
clear

x1=[120 140 190 130 155 175 125 145 180 150];
x2=[100 110  90 150 210 150 250 270 300 250];
y= [102 100 120  77  46  93  26  69  65  85]';
x=[ones(10,1) x1' x2'];
x_O=[x1' x2'];

[b,bint,r,rint,stats]=regress(y,x_O)
figure;rcoplot(r,rint);