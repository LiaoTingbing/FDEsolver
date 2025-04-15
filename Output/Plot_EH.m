
clc ;
close all;
clear ;

Ex_real = load("Ex_real.txt");
Ex_imag = load("Ex_imag.txt");

Ey_real = load("Ey_real.txt");
Ey_imag = load("Ey_imag.txt");

Ez_real = load("Ez_real.txt");
Ez_imag = load("Ez_imag.txt");

Hx_real = load("Hx_real.txt");
Hx_imag = load("Hx_imag.txt");

Hy_real = load("Hy_real.txt");
Hy_imag = load("Hy_imag.txt");

Hz_real = load("Hz_real.txt");
Hz_imag = load("Hz_imag.txt");

x = load("x.txt");
y = load("y.txt");

neff_real = load("neff_real.txt");
neff_imag = load("neff_imag.txt");

Ex = Ex_real + 1i * Ex_imag;
Ey = Ey_real + 1i * Ey_imag;
Ez = Ez_real + 1i * Ez_imag;

Hx = Hx_real + 1i * Hx_imag;
Hy = Hy_real + 1i * Hy_imag;
Hz = Hz_real + 1i * Hz_imag;


neff = neff_real + 1i * neff_imag;

Ex = reshape(Ex , length(x) , length(y) , []);
Ey = reshape(Ey , length(x) , length(y) , []);
Ez = reshape(Ez , length(x) , length(y) , []);

Hx = reshape(Hx , length(x) , length(y) , []);
Hy = reshape(Hy , length(x) , length(y) , []);
Hz = reshape(Hz , length(x) , length(y) , []);


Eamp  =  sqrt(abs(Ex).^2 + abs(Ey).^2 + abs(Ey).^2 );
Hamp  =  sqrt(abs(Hx).^2 + abs(Hy).^2 + abs(Hy).^2 );
%% 
 
[X,Y]  =ndgrid(x,y);
mode = 3;

% figure
% pcolor(X,Y,(abs(Ex(:,:,mode))));subtitle("Ex");colorbar;
% axis equal;shading interp;colormap jet;
% xlabel("x (um)") ; ylabel("y (um)")
% 
% 
% figure
% pcolor(X,Y,abs(Ey(:,:,mode))); subtitle("Ey");colorbar;
% axis equal;shading interp;colormap jet
% xlabel("x (um)") ; ylabel("y (um)")
% 
% figure
% pcolor(X,Y,abs(Ez(:,:,mode)));subtitle("Ez");colorbar;
% axis equal;shading interp;colormap jet
%  xlabel("x (um)") ; ylabel("y (um)")
% 
% figure
% pcolor(X,Y,abs(Hx(:,:,mode)));subtitle("Hx");colorbar;
% axis equal;shading interp;colormap jet
% xlabel("x (um)") ; ylabel("y (um)")
% 
% figure
% pcolor(X,Y,abs(Hy(:,:,mode)));subtitle("Hy");colorbar;
% axis equal;shading interp;colormap jet
% xlabel("x (um)") ; ylabel("y (um)")
% 
% figure
% pcolor(X,Y,abs(Hz(:,:,mode)));subtitle("Hz");colorbar;
% axis equal;shading interp;colormap jet
% xlabel("x (um)") ; ylabel("y (um)")

%%

% subplot(1,2,1)
figure
pcolor(X,Y,abs(Eamp(:,:,mode)));subtitle("E intensity");colorbar;
axis equal;shading interp;colormap jet
xlabel("x (um)") ; ylabel("y (um)")

figure
% subplot(1,2,2)
pcolor(X,Y,abs(Hamp(:,:,mode)));subtitle("H intensity");colorbar;
axis equal;shading interp;colormap jet
xlabel("x (um)") ; ylabel("y (um)")
 


 