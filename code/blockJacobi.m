function [x,err] = Bjacobi(A,b,x,xstar,MaxIt,tol)
%
% Solve Ax = b using Block Jacobi iteration
%
% with initial guess x
%
b = b(:); x = x(:); err = zeros(MaxIt,1);
DB = [A(1:2,1:2) zeros(2,2); zeros(2,2) A(3:4,3:4)]; % Block diag of A
disp(sprintf(’x = ’));
for k = 1:MaxIt
r = b - A*x;
if (norm(r)<tol*norm(b)), break; end
x = x + DB\r;
err(k) = norm(x-xstar);
disp(sprintf(’%12.4e ’, x));
end