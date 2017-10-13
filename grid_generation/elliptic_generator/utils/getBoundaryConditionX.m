
% @brief: returns the value at the boundary condition given by the position i, j for the ''x'' component
% @param {int} i : current ith position ( where stencil is evaluated )
% @param {int} j : current jth position ( where stencil is evaluated )
% @param {int} in : current ith position of the boundary given the stencil
% @param {int} jn : current jth position of the boundary given the stencil
% @param {int} gSz_i : size of the grid in the i component
% @param {int} gSz_j : size of the grid in the j component
% @return {float} valBoundary : value at this boundary point
function valBoundary = getBoundaryConditionX( xg, yg, i, j, in, jn, gSz_i, gSz_j )

    di = in - i;
    dj = jn - j;

    if di == 1 && dj == 1

        valBoundary =  -0.5 * c_beta( xg, yg, i, j ) * xg( in + 1, jn + 1 );

    elseif di == 1 && dj == 0

        valBoundary =        c_alpha( xg, yg, i, j ) * xg( in + 1, jn + 1 );

    elseif di == 1 && dj == -1

        valBoundary =   0.5 * c_beta( xg, yg, i, j ) * xg( in + 1, jn + 1 );

    elseif di == 0 && dj == -1

        valBoundary =        c_gamma( xg, yg, i, j ) * xg( in + 1, jn + 1 );

    elseif di == -1 && dj == -1

        valBoundary = -0.5 * c_beta( xg, yg, i , j ) * xg( in + 1, jn + 1 );

    elseif di == -1 && dj == 0

        valBoundary =        c_alpha( xg, yg, i, j ) * xg( in + 1, jn + 1 );

    elseif di == -1 && dj == 1

        valBoundary =   0.5 * c_beta( xg, yg, i, j ) * xg( in + 1, jn + 1 );

    elseif di == 0 && dj == 1

        valBoundary =        c_gamma( xg, yg, i, j ) * xg( in + 1, jn + 1 );

    else 
        disp( 'shouldnt get here ' );
    end


end