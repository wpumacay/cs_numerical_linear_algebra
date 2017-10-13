function [x,err] = gauss_seidel_solver(A,b,x,xstar,MaxIt,tol)
	%
	% Solve Ax = b using Gauss-Siedel iteration
	%
	% 
	%
	b = b(:); x = x(:); err = zeros(MaxIt,1);
	DL = tril(A); U = -triu(A,1);
	%disp(sprintf(’x = ’));
	for k = 1:MaxIt
	r = b - A*x;
	if (norm(r)<tol*norm(b)), break; end
	x = x + DL\r;
	err(k) = norm(x-xstar);
	%disp(sprintf(’%12.4e ’, x));
	end
	disp(sprintf('Number of iterations: %d ', k));
end
