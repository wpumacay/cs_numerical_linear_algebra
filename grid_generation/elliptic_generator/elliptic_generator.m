
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
   	new_xg = basicSysSolver( Ax, bx );
   	new_yg = basicSysSolver( Ay, by );

    xg( 2 : end - 1, 2 : end - 1 ) = reshape( new_xg, grid_sz_i, grid_sz_j );
    yg( 2 : end - 1, 2 : end - 1 ) = reshape( new_yg, grid_sz_i, grid_sz_j );

end