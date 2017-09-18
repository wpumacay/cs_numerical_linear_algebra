
% @brief: return the discretized system matrix and vector ( A, b ) for both x and y coordinates
% @param {matrix} xg : x-grid
% @param {matrix} yg : y-grid
% @param {int} tSz_i : size of the entire grid in the i component ( including boundaries )
% @param {int} tSz_j : size of the entire grid in the j component ( including boundaries )
% @param {int} gSz_i : size of the grid in the i component ( just part to apply stencil on )
% @param {int} gSz_j : size of the grid in the j component ( just part to apply stencil on )
% @return {[matrix, matrix]} [A, b] : system's matrix and vector
function [Ax, Ay, bx, by] = generateSysMatrix( xg, yg, tSz_i, tSz_j, gSz_i, gSz_j )
    dim = gSz_i * gSz_j;
    Ax = zeros( dim, dim );
    Ay = zeros( dim, dim );

    bx = zeros( dim, 1 );
    by = zeros( dim, 1 );
    
    for i = 1 : gSz_i
        for j = 1 : gSz_j
            [Ax, bx] = applyStencilX( Ax, bx, xg, yg, i, j, gSz_i, gSz_j );
            [Ay, by] = applyStencilY( Ay, by, xg, yg, i, j, gSz_i, gSz_j );
        end
    end

end