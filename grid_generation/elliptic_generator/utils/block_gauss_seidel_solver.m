function [x,err] = block_gauss_seidel_solver(A,b,x,xstar,MaxIt,tol)
	%
	% Solve Ax = b using Block Gauss-Siedel iteration
	%
	% with initial guess x
	% 
	b = b(:); 
    x = x(:); 
    err = zeros( MaxIt, 1 );
	% DL = [A(1:2,1:2) zeros(2,2); A(3:4,1:2) A(3:4,3:4)]; U = A - DL;
	%DL = [A(1:5,1:5) zeros(5,4); A(6:9,1:5) A(6:9,6:9)]; U = A - DL;
	m = round( length( b ) / 2 );
	n = length( b ) - round( length( b ) / 2 );
	DL = [A(1:m,1:m) zeros(m,n); A(m+1:m+n,1:m) A(m+1:m+n,m+1:m+n)]; U = A - DL;
	%disp( 'x = ' );
	for k = 1:MaxIt
	  r = b - A*x;
	  if (norm(r)<tol*norm(b)), break; end
	  x = x + DL\r;
	  err(k) = norm(x-xstar);
	  %disp(sprintf('%12.4e ', x));
	end
	fprintf( 'Number of iterations: %d ', k );
end
