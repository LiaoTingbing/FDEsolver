
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
modeselect = 2;       % 基模1
componentselect = 2 ;     % 1为Ex % 2为Ey
fieldcomponent = {"Ex" , "Ey"};
start = (componentselect-1)*nx*ny+1;
for i=modeselect
pcolor(reshape(field.real(start:start+nx*ny-1,i),nx,[]))
shading interp
title( "模式" + num2str(modeselect) + "   "+ fieldcomponent{componentselect}+ "   有效折射率   "+num2str(neff.real(i)))
colormap jet
colorbar
pause(1)
axis equal
saveas(gcf,"../doc/images/fiber/mode"+num2str(modeselect)+fieldcomponent{componentselect}+".png")
end