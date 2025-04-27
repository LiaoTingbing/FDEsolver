

function [s]=dydyfuncVec(p,q,r,nx,ny,dx,dy)

p =1./p;
ng = nx * ny;
z  = zeros(nx,1);
a =  p.*[z;r(1:ng-nx)]./([z;q(1:ng-nx)] + q);
b = - p.*r.*( 1./([z;q(1:ng-nx)] + q) + 1./(q +[q(nx+1:ng);z]) );
c =  p.*[r(nx+1:ng);z]./(q +[q(nx+1:ng);z]);

s =  transpose( spdiags([c,b,a]*2/dy/dy,[-nx 0 nx],ng,ng));
end


%%
% 函数测试
% clear
% ss = dydyfunc( (1:16)',(1:16)',(1:16)',4,4,1,1);
% 
% 
% spy(ss)