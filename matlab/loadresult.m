
clear;

filename = "out.h5" ; 

x = h5read(filename , '/x');
y = h5read(filename , '/y');
nx = length(x); 
ny = length(y);


field  = h5read(filename , '/field');
neff  = h5read(filename , '/neff');
 

%%
close
for i=1:40
pcolor(reshape(field.real(1:nx*ny,i),nx,[]))
shading interp
colormap jet
colorbar
pause(1)
axis equal
end