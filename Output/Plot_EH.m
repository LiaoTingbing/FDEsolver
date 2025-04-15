
clc ;
close all;
clear ;

model = "waveguide";
% model = "fiber";
inPath = model + "\";

%%  LOAD DATA

Ex_real = load(inPath+"Ex_real.txt");
Ex_imag = load(inPath+"Ex_imag.txt");

Ey_real = load(inPath+"Ey_real.txt");
Ey_imag = load(inPath+"Ey_imag.txt");

Ez_real = load(inPath+"Ez_real.txt");
Ez_imag = load(inPath+"Ez_imag.txt");

Hx_real = load(inPath+"Hx_real.txt");
Hx_imag = load(inPath+"Hx_imag.txt");

Hy_real = load(inPath+"Hy_real.txt");
Hy_imag = load(inPath+"Hy_imag.txt");

Hz_real = load(inPath+"Hz_real.txt");
Hz_imag = load(inPath+"Hz_imag.txt");

x = load(inPath+"x.txt");
y = load(inPath+"y.txt");

neff_real = load(inPath+"neff_real.txt");
neff_imag = load(inPath+"neff_imag.txt");

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



for i =progress(1:5) 
    % mode = 3;
    mode = i;
    outPath = inPath +"mode" + num2str(mode) + "\";

    %%
    system("mkdir " + outPath  );
    [X,Y]  =ndgrid(x,y);
    sf = { "Ex" , "Ey" , "Ez" , "Hx" , "Hy" , "Hz" };
    for i=1:6
        figure
        f = eval( " abs(" + sf{i} + "(:,:,mode)) " );
        pcolor(X,Y,f);subtitle(sf{i});colorbar;
        title("有效折射率：" + num2str(neff(mode)))
        axis equal;shading interp;colormap jet
        xlabel("x (um)") ; ylabel("y (um)")
        saveas(gcf  ,outPath+ "mode" + num2str(mode)+"_"+ sf{i}+"_amplitude.png");
    end

    %% EH intensity

    figure;
    pcolor(X,Y,abs(Eamp(:,:,mode)));subtitle("E intensity");colorbar;
    axis equal;shading interp;colormap jet
    title("有效折射率：" + num2str(neff(mode)))
    xlabel("x (um)") ; ylabel("y (um)")
    saveas( gcf , outPath+ "mode" + num2str(mode)+"_E_intensity.png")

    figure
    pcolor(X,Y,abs(Hamp(:,:,mode)));subtitle("H intensity");colorbar;
    axis equal;shading interp;colormap jet
    title("有效折射率：" + num2str(neff(mode)))
    xlabel("x (um)") ; ylabel("y (um)")
    saveas( gcf , outPath+ "mode" + num2str(mode)+"_H_intensity.png")

end

%% 
close all