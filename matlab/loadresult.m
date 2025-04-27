
clear;
result = readmatrix("..\text.txt");

pcolor(abs(reshape(result(1:101*101,15),101,101)));
shading interp
colormap jet