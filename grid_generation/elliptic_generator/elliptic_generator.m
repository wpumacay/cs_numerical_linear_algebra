
% @brief: returns a refined version of the x-y grid
% @param {matrix} xg : x-grid
% @param {matrix} yg : y-grid
% @return {[matrix, matrix]} [xg, yg] : refined x-y grid
function [xg, yg] = elliptic_generator( xg, yg )

	sz_i = size( xg, 1 );
	sz_j = size( xg, 2 );

	grid_sz_i = sz_i - 2;
	grid_sz_j = sz_j - 2;

    [Ax, Ay, bx, by] = generateSysMatrix( xg, yg, sz_i, sz_j, ...
                                          grid_sz_i, grid_sz_j );

    % Replace this basic solver by an iterative method
    xg_0 = xg( 2 : end - 1, 2 : end - 1 ); 
    yg_0 = yg( 2 : end - 1, 2 : end - 1 );
    xg_0 = reshape( xg_0, size( xg_0, 1 ) * size( xg_0, 2 ), 1 );
    yg_0 = reshape( yg_0, size( yg_0, 1 ) * size( yg_0, 2 ), 1 );

    N = size( Ax, 1 );

    xstar = Ax \ bx; 
    ystar = Ay \ by;

    maxIt = 2000; 
    tol = 1e-6;

   	[new_xg, ~] = block_gauss_seidel_solver( Ax, bx, xg_0, xstar, maxIt, tol );
   	[new_yg, ~] = block_gauss_seidel_solver( Ay, by, yg_0, ystar, maxIt, tol );

    %[new_xg, ~] = gauss_seidel_solver( Ax, bx, xg_0, xstar, maxIt, tol );
    %[new_yg, ~] = gauss_seidel_solver( Ay, by, yg_0, ystar, maxIt, tol );

    %[new_xg, t1] = basicSysSolver( Ax, bx );
    %[new_yg, t2] = basicSysSolver( Ay, by );

    xg( 2 : end - 1, 2 : end - 1 ) = reshape( new_xg, grid_sz_i, grid_sz_j );
    yg( 2 : end - 1, 2 : end - 1 ) = reshape( new_yg, grid_sz_i, grid_sz_j );

end