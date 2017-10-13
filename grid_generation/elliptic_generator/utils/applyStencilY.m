
% @brief: return the discretized system matrix and vector ( A, b ) to be resolved to refine the grid
% @param {matrix} A : system matrix of the discretization
% @param {matrix} b : system vector of the discretization
% @param {matrix} xg : x-grid
% @param {matrix} yg : y-grid
% @param {int} i : current ith position ( where stencil is evaluated )
% @param {int} j : current jth position ( where stencil is evaluated )
% @param {int} gSz_i : size of the grid in the i component
% @param {int} gSz_j : size of the grid in the j component
% @return {[matrix, matrix]} [A, b] : system's matrix and vector
function [A, b] = applyStencilY( A, b, xg, yg, i, j, gSz_i, gSz_j )

    % check if is corner or border or interior
    
    s_row = getIndxInMatrix( i, j, gSz_i, gSz_j );
    
    if i == 1 % in bottom border
        
        if j == 1 % bottom left corner
            
            center_col   = s_row;
            right_col    = getIndxInMatrix( i, j + 1, gSz_i, gSz_j );
            right_up_col = getIndxInMatrix( i + 1, j + 1, gSz_i, gSz_j );
            up_col       = getIndxInMatrix( i + 1, j, gSz_i, gSz_j );
            
            A( s_row, center_col )      = -2 * c_alpha( xg, yg, i, j ) - ...
                                           2 * c_gamma( xg, yg, i, j );
            A( s_row, right_col )       = c_gamma( xg, yg, i, j );
            A( s_row, right_up_col )    = -0.5 * c_beta( xg, yg, i, j );
            A( s_row, up_col )          = c_alpha( xg, yg, i, j );
            
            b( s_row, 1 ) = -getBoundaryConditionY( xg, yg, i, j, i + 1, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i - 1, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i - 1, j, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i - 1, j + 1, gSz_i, gSz_j );
            
        elseif j == gSz_j % bottom right corner
            
            center_col      = s_row;
            left_up_col     = getIndxInMatrix( i + 1, j - 1, gSz_i, gSz_j );
            left_col        = getIndxInMatrix( i, j - 1, gSz_i, gSz_j );
            up_col          = getIndxInMatrix( i + 1, j, gSz_i, gSz_j );

            A( s_row, center_col )      = -2 * c_alpha( xg, yg, i, j ) - ...
                                           2 * c_gamma( xg, yg, i, j );
            A( s_row, left_up_col )     = 0.5 * c_beta( xg, yg, i, j );
            A( s_row, left_col )        = c_gamma( xg, yg, i, j );
            A( s_row, up_col )          = c_alpha( xg, yg, i, j );
            
            b( s_row, 1 ) = -getBoundaryConditionY( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i    , j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i - 1, j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i - 1, j    , gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i - 1, j - 1, gSz_i, gSz_j );
            
        else % a point in the bottom border, but not corner

            center_col      = s_row;
            left_up_col     = getIndxInMatrix( i + 1, j - 1, gSz_i, gSz_j );
            left_col        = getIndxInMatrix( i, j - 1, gSz_i, gSz_j );
            right_col       = getIndxInMatrix( i, j + 1, gSz_i, gSz_j );
            right_up_col    = getIndxInMatrix( i + 1, j + 1, gSz_i, gSz_j );
            up_col          = getIndxInMatrix( i + 1, j, gSz_i, gSz_j );

            A( s_row, center_col )      = -2 * c_alpha( xg, yg, i, j ) - ...
                                           2 * c_gamma( xg, yg, i, j );
            A( s_row, left_up_col )     = 0.5 * c_beta( xg, yg, i, j );
            A( s_row, left_col )        = c_gamma( xg, yg, i, j );
            A( s_row, right_col )       = c_gamma( xg, yg, i, j );
            A( s_row, right_up_col )    = -0.5 * c_beta( xg, yg, i, j );
            A( s_row, up_col )          = c_alpha( xg, yg, i, j );
            
            b( s_row, 1 ) = -getBoundaryConditionY( xg, yg, i, j, i - 1, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i - 1, j    , gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i - 1, j + 1, gSz_i, gSz_j );
            
        end
        
    elseif i == gSz_i
        
        if j == 1 % top left corner

            center_col      = s_row;
            down_col        = getIndxInMatrix( i - 1, j, gSz_i, gSz_j );
            right_down_col  = getIndxInMatrix( i - 1, j + 1, gSz_i, gSz_j );
            right_col       = getIndxInMatrix( i, j + 1, gSz_i, gSz_j );

            A( s_row, center_col )      = -2 * c_alpha( xg, yg, i, j ) - ...
                                           2 * c_gamma( xg, yg, i, j );
            A( s_row, down_col )        = c_alpha( xg, yg, i, j );
            A( s_row, right_down_col )  = 0.5 * c_beta( xg, yg, i, j );
            A( s_row, right_col )       = c_gamma( xg, yg, i, j );
            
            b( s_row, 1 ) = -getBoundaryConditionY( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i + 1, j    , gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i + 1, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i    , j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i - 1, j - 1, gSz_i, gSz_j );
            
        elseif j == gSz_j % topright corner

            center_col      = s_row;
            left_col        = getIndxInMatrix( i, j - 1, gSz_i, gSz_j );
            left_down_col   = getIndxInMatrix( i - 1, j - 1, gSz_i, gSz_j );
            down_col        = getIndxInMatrix( i - 1, j, gSz_i, gSz_j );

            A( s_row, center_col )      = -2 * c_alpha( xg, yg, i, j ) - ...
                                           2 * c_gamma( xg, yg, i, j );
            A( s_row, left_col )        = c_gamma( xg, yg, i, j );
            A( s_row, left_down_col )   = -0.5 * c_beta( xg, yg, i, j );
            A( s_row, down_col )        = c_alpha( xg, yg, i, j );
            
            b( s_row, 1 ) = -getBoundaryConditionY( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i + 1, j    , gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i    , j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i - 1, j + 1, gSz_i, gSz_j );
            
        else % a point in the top border, but not corner

            center_col      = s_row;
            left_col        = getIndxInMatrix( i, j - 1, gSz_i, gSz_j );
            left_down_col   = getIndxInMatrix( i - 1, j - 1, gSz_i, gSz_j );
            down_col        = getIndxInMatrix( i - 1, j, gSz_i, gSz_j );
            right_down_col  = getIndxInMatrix( i - 1, j + 1, gSz_i, gSz_j );
            right_col       = getIndxInMatrix( i, j + 1, gSz_i, gSz_j );

            A( s_row, center_col )      = -2 * c_alpha( xg, yg, i, j ) - ...
                                           2 * c_gamma( xg, yg, i, j );
            A( s_row, left_col )        = c_gamma( xg, yg, i, j );
            A( s_row, left_down_col )   = -0.5 * c_beta( xg, yg, i, j );
            A( s_row, down_col )        = c_alpha( xg, yg, i, j );
            A( s_row, right_down_col )  = 0.5 * c_beta( xg, yg, i, j );
            A( s_row, right_col )       = c_gamma( xg, yg, i, j );
            
            b( s_row, 1 ) = -getBoundaryConditionY( xg, yg, i, j, i + 1, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i + 1, j    , gSz_i, gSz_j ) - ...
                             getBoundaryConditionY( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j );            
            
        end
        
    elseif j == 1
        
        % left border, but not corner

        center_col      = s_row;
        down_col        = getIndxInMatrix( i - 1, j, gSz_i, gSz_j );
        right_down_col  = getIndxInMatrix( i - 1, j + 1, gSz_i, gSz_j );
        right_col       = getIndxInMatrix( i, j + 1, gSz_i, gSz_j );
        right_up_col    = getIndxInMatrix( i + 1, j + 1, gSz_i, gSz_j );
        up_col          = getIndxInMatrix( i + 1, j, gSz_i, gSz_j );
        
        A( s_row, center_col )      = -2 * c_alpha( xg, yg, i, j ) - ...
                                       2 * c_gamma( xg, yg, i, j );
        A( s_row, down_col )        = c_alpha( xg, yg, i, j );
        A( s_row, right_down_col )  = 0.5 * c_beta( xg, yg, i, j );
        A( s_row, right_col )       = c_gamma( xg, yg, i, j );
        A( s_row, right_up_col )    = -0.5 * c_beta( xg, yg, i, j );
        A( s_row, up_col )          = c_alpha( xg, yg, i, j );
        
        b( s_row, 1 ) = -getBoundaryConditionY( xg, yg, i, j, i + 1, j - 1, gSz_i, gSz_j ) - ...
                         getBoundaryConditionY( xg, yg, i, j, i    , j - 1, gSz_i, gSz_j ) - ...
                         getBoundaryConditionY( xg, yg, i, j, i - 1, j - 1, gSz_i, gSz_j );            
        
    elseif j == gSz_j
        
        % right border, but not corner

        center_col      = s_row;
        left_up_col     = getIndxInMatrix( i + 1, j - 1, gSz_i, gSz_j );
        left_col        = getIndxInMatrix( i, j - 1, gSz_i, gSz_j );
        left_down_col   = getIndxInMatrix( i - 1, j - 1, gSz_i, gSz_j );
        down_col        = getIndxInMatrix( i - 1, j, gSz_i, gSz_j );
        up_col          = getIndxInMatrix( i + 1, j, gSz_i, gSz_j );
        
        A( s_row, center_col )      = -2 * c_alpha( xg, yg, i, j ) - ...
                                       2 * c_gamma( xg, yg, i, j );
        A( s_row, left_up_col )     = 0.5 * c_beta( xg, yg, i, j );
        A( s_row, left_col )        = c_gamma( xg, yg, i, j );
        A( s_row, left_down_col )   = -0.5 * c_beta( xg, yg, i, j );
        A( s_row, down_col )        = c_alpha( xg, yg, i, j );
        A( s_row, up_col )          = c_alpha( xg, yg, i, j );
        
        b( s_row, 1 ) = -getBoundaryConditionY( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j ) - ...
                         getBoundaryConditionY( xg, yg, i, j, i    , j + 1, gSz_i, gSz_j ) - ...
                         getBoundaryConditionY( xg, yg, i, j, i - 1, j + 1, gSz_i, gSz_j );            
    else
        
        % interior point
        center_col      = s_row;
        left_up_col     = getIndxInMatrix( i + 1, j - 1, gSz_i, gSz_j );
        left_col        = getIndxInMatrix( i, j - 1, gSz_i, gSz_j );
        left_down_col   = getIndxInMatrix( i - 1, j - 1, gSz_i, gSz_j );
        down_col        = getIndxInMatrix( i - 1, j, gSz_i, gSz_j );
        right_down_col  = getIndxInMatrix( i - 1, j + 1, gSz_i, gSz_j );
        right_col       = getIndxInMatrix( i, j + 1, gSz_i, gSz_j );
        right_up_col    = getIndxInMatrix( i + 1, j + 1, gSz_i, gSz_j );
        up_col          = getIndxInMatrix( i + 1, j, gSz_i, gSz_j );
        
        A( s_row, center_col )      = -2 * c_alpha( xg, yg, i, j ) - ...
                                       2 * c_gamma( xg, yg, i, j );
        A( s_row, left_up_col )     = 0.5 * c_beta( xg, yg, i, j );
        A( s_row, left_col )        = c_gamma( xg, yg, i, j );
        A( s_row, left_down_col )   = -0.5 * c_beta( xg, yg, i, j );
        A( s_row, down_col )        = c_alpha( xg, yg, i, j );
        A( s_row, right_down_col )  = 0.5 * c_beta( xg, yg, i, j );
        A( s_row, right_col )       = c_gamma( xg, yg, i, j );
        A( s_row, right_up_col )    = -0.5 * c_beta( xg, yg, i, j );
        A( s_row, up_col )          = c_alpha( xg, yg, i, j );
        
    end
    
end