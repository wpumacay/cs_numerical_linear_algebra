function [x,N]= BJ(A,b,x0,d,eps,M) %块雅克比迭代法求线性方程组Ax=b的解

if nargin==4

eps= 1.0e-6;

M = 10000;

elseif nargin<4

error

return

elseif nargin ==5

M = 10000; %参数的默认值

end

NS = size(A);

n = NS(1,1);

if(sum(d) ~= n)

disp('分块错误！');

return;

end

bnum = length(d);

bs = ones(bnum,1);

for i=1:(bnum-1)

bs(i+1,1)=sum(d(1:i))+1;