
% @brief: return the alpha coefficient of the discretization of the ...
% pde physical coordinates respecto to curvilinear coordinates
% @param {matrix} xg : x-grid
% @param {matrix} yg : y-grid
% @param {int} i : current ith position to take the coefficient ( where stencil is evaluated )
% @param {int} j : current jth position to take the coefficient ( where stencil is evaluated )
% @return {float} g : alpha coefficient at the given position
function g = c_alpha( xg, yg, i, j )
    g = 0.25 * ( ( xg( i + 1, j + 2 ) - xg( i + 1, j ) ) ^ 2 + ...
                 ( yg( i + 1, j + 2 ) - yg( i + 1, j ) ) ^ 2 );
end