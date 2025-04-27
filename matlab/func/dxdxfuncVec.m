
% p q r 系数
function [s]=dxdxfuncVec(p,q,r,nx,ny,dx,dy)

p =1./p;

ng = nx * ny;
a =  p.*[0;r(1:ng-1)]./([0;q(1:ng-1)]+q);
b = - p.*r.*( 1./([0;q(1:ng-1)]+q) + 1./(q + [q(2:ng);0]) ) ;
c =  p.*[r(2:ng);0]./( q + [q(2:ng);0]) ;

a(1:nx:ng)=0;
c(nx:nx:ng)=0;

s = transpose( spdiags([c b a]*2/dx/dx , [-1 0 1] , ng ,ng ) );

end
%%
% 函数测试
% ss = dxdxfuncVec( (1:16)',(1:16)',(1:16)',4,4,1,1)
% 
% spy(ss)