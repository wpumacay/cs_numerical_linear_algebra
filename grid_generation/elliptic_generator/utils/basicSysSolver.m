
% @brief : solves Ax=b using matlab's '\' operator
function [x,t] = basicSysSolver( A, b )
	tic;
	x = A \ b;
	t = toc;
end