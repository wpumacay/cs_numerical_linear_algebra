function [x,err] = Bgs(A,b,x,xstar,MaxIt,tol)
%
% Solve Ax = b using Block Gauss-Siedel iteration
%
b = b(:); x = x(:); err = zeros(MaxIt,1);
m = round(size(b)/2);
n = size(b) - round(size(b)/2);
DL = [A(1:m,1:m) zeros(m,n); A(m+1:m+n,1:m) A(m+1:m+n,m+1:m+n)]; U = A - DL;
disp(sprintf('x = '));
for k = 1:MaxIt
  r = b - A*x;
  if (norm(r)<tol*norm(b)), break; end
  x = x + DL\r;
  err(k) = norm(x-xstar);
end
disp(sprintf('Number of iterations: %d ', k));
end

