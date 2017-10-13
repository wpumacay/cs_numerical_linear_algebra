
% @brief: returns the indx of the current i,j position in the global matrix that represents the discretized system
% @param {int} i : current ith position ( where stencil is evaluated )
% @param {int} j : current jth position ( where stencil is evaluated )
% @param {int} gSz_i : size of the grid in the i component
% @return {indx} indx : index in the system matrix
function indx = getIndxInMatrix( i, j, gSz_i, ~ )
    indx = i + ( j - 1 ) * gSz_i;
end