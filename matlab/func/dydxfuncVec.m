

function [s] = dydxfuncVec(p,q,r,nx,ny,dx,dy)

p =1./p;
ng = nx * ny ;

a =   p.*[zeros(nx+1,1);r(1:ng-(nx+1))].*[zeros(nx,1);1./q(1:ng-nx)];
b = - p.*[zeros(nx-1,1);r(1:ng-(nx-1))].*[zeros(nx,1);1./q(1:ng-nx)];
c = - p.*[r(nx:ng);zeros(nx-1,1)].*[1./q(nx+1:ng);zeros(nx,1)];
d =   p.*[r(nx+2:ng);zeros(nx+1,1)].*[1./q(nx+1:ng);zeros(nx,1)];

a(1:nx:ng)=0;
b(nx:nx:ng)=0;
c(1:nx:ng)=0;
d(nx:nx:ng)=0;

s = transpose( spdiags( [d c b a]/4/dx/dy ,...
    [-nx-1 , -nx+1 , nx-1 , nx+1] , ng,ng) );
 
end



%%
% % 函数测试
% ss =dydxfuncVec1( (1:16)',(1:16)',(1:16)',4,4,1,1)
%
% spy(ss)