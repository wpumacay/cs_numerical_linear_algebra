
% @brief: return the gamma coefficient of the discretization of the ...
% pde physical coordinates respecto to curvilinear coordinates
% @param {matrix} xg : x-grid
% @param {matrix} yg : y-grid
% @param {int} i : current ith position to take the coefficient ( where stencil is evaluated )
% @param {int} j : current jth position to take the coefficient ( where stencil is evaluated )
% @return {float} g : gamma coefficient at the given position
function g = c_beta( xg, yg, i, j )
    g = 0.25 * ( ( xg( i + 2, j + 1 ) - xg( i, j + 1 ) ) ^ 2 + ...
                 ( yg( i + 2, j + 1 ) - yg( i, j + 1 ) ) ^ 2 );
end