
%基于BPM
clear ;

eps0 =8.85419e-12; 

filename = "../lumerical/lumerical.h5";

x = h5read(filename , "/x");
y = h5read(filename , "/y");

dx_ = x(2)-x(1);
dy_ = y(2)-y(1);
nx_ = length(x);
ny_ = length(y);
nt_ = nx_ * ny_ ;

lambda = h5read(filename , "/lambda");
k0_ = 2*pi/lambda;

indexX = h5read(filename , "/indexX");
indexY = h5read(filename , "/indexY");
indexZ = h5read(filename , "/indexZ");
indexXY = h5read(filename , "/indexXY");
indexYX = h5read(filename , "/indexYX");


epsZ = indexZ(:);
epsX = indexX(:);
epsY = indexY(:);
epsXY = indexXY(:);
epsYX = indexYX(:);

%% PML 参数 10 层
c0 = 3e8;
Impedance = 373;
w0 = k0_ * c0;
R0 = 1e-16 ;   % 反射率
SigmaXMax = -( 3 + 1  ) * log(R0) / 2 / Impedance / (10*dx_)     ;
SigmaYMax = -( 3 + 1  ) * log(R0) / 2 / Impedance / (10*dy_)     ;

sx = ones( nx_ , ny_ ) ;
% sx(1:11 , : ) = repmat(  1 + (( 10:-1:0)/10 ).' .^3 * SigmaXMax  / 1j / w0  /eps0  , 1 , ny_) ;
% sx(end:-1:end-10 , : ) = repmat(  1 + (( 10:-1:0)/10 ).' .^3 * SigmaXMax  / 1j / w0  /eps0  , 1 , ny_) ;

sy  = ones( nx_ , ny_ ) ;
% sy(: , 1:11 ) = repmat(  1 + (( 10:-1:0)/10 )  .^3 * SigmaXMax  / 1j / w0  /eps0  , nx_ , 1) ;
% sy(:  , end:-1:end-10 ) = repmat(  1 + (( 10:-1:0)/10 )  .^3 * SigmaXMax  / 1j / w0  /eps0  , nx_ , 1) ;
% 
% sxvec =  sx(:);
% syvec =  sy(:);

sx = sx(:);
sy = sy(:);
isx = 1./sx;
isy = 1./sy;

%%
Pxx_ = dxdxfuncVec(isx,sx.*epsZ, epsX, nx_, ny_, dx_, dy_)...
+ dydyfuncVec(isy, sy,ones(nt_,1), nx_, ny_, dx_, dy_)...
+ spdiags(k0_ * k0_ * epsX, 0, nt_, nt_)...
+ dxdyfuncVec(isx, sy.*epsZ, epsYX, nx_, ny_, dx_, dy_);

Pxy_ = dxdyfuncVec(isx,sy.* epsZ, epsY, nx_, ny_, dx_, dy_)...
- dxdyfuncVec(isx,sy, ones(nt_,1), nx_, ny_, dx_, dy_)...
+ spdiags(k0_ * k0_ * epsXY, 0, nt_, nt_)...
+ dxdxfuncVec(isx, sx.*epsZ, epsXY, nx_, ny_, dx_, dy_);

Pyx_ = dydxfuncVec(isy, sx.*epsZ, epsX, nx_, ny_, dx_, dy_)...
- dydxfuncVec(isy, sx, ones(nt_,1), nx_, ny_, dx_, dy_)...
+ spdiags(k0_ * k0_ * epsYX, 0, nt_, nt_)...
+ dydyfuncVec(isy, sy.*epsZ, epsYX, nx_, ny_, dx_, dy_);

Pyy_ = dydyfuncVec(isy, sy.*epsZ, epsY, nx_, ny_, dx_, dy_)...
+ dxdxfuncVec(isx, sx, ones(nt_,1), nx_, ny_, dx_, dy_)...
+ spdiags(k0_ * k0_ * epsY, 0, nt_, nt_)...
+ dydxfuncVec(isy, sx.*epsZ, epsXY, nx_, ny_, dx_, dy_);


P = [Pxx_ , Pxy_ ; Pyx_ , Pyy_];

[e ,v] = eigs(P,20,3.4*3.4*k0_*k0_);

neff = sqrt(diag(v))/k0_;

ex = reshape(e(1:nx_*ny_,:),nx_,ny_,[]);
ey = reshape(e( nx_*ny_+1:end,:),nx_,ny_,[]);
%%
% close
% for i=1:10
% pcolor(abs(ex(:,:,i)))
% title(num2str(neff(i)))
% axis equal
% shading interp
% colormap jet
% colorbar
% pause(1)
% end
%%
