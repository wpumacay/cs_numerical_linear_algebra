function [x,err] = Bgs(A,b,x,xstar,MaxIt,tol)
%
% Solve Ax = b using Block Gauss-Siedel iteration
%
% with initial guess x
%
b = b(:); x = x(:); err = zeros(MaxIt,1);
DL = [A(1:2,1:2) zeros(2,2); A(3:4,1:2) A(3:4,3:4)]; U = A - DL;
disp(sprintf('x = '));
for k = 1:MaxIt
r = b - A*x;
if (norm(r)<tol*norm(b)), break; end
x = x + DL\r;
err(k) = norm(x-xstar);
disp(sprintf('%12.4e ', x));
end
disp(sprintf('Number of iterations: %d ', k));