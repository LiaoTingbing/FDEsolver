
clear;

filename = "out.h5" ; 

x = h5read(filename , '/x');
y = h5read(filename , '/y');
nx = length(x); 
ny = length(y);

fieldAbs = h5read(filename , '/fieldAbs');

pcolor(reshape(fieldAbs(1:nx*ny,5),nx,[]))
shading interp
colormap jet