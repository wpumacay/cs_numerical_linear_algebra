disp( 'running algebraic generator' );
alg_geometry_loader;
disp( 'finished running algebraic generator' );

pause

sz_i = size( xx_g, 1 );
sz_j = size( xx_g, 2 );

grid_sz_i = sz_i - 2;
grid_sz_j = sz_j - 2;


Niter = 10;

for k = 1 : Niter
   
    [Ax, Ay, bx, by] = generateSysMatrix( xx_g, yy_g, sz_i, ...
                                          sz_j, grid_sz_i, grid_sz_j );
%     disp( 'Ax' )
%     disp( Ax )
%     disp( 'Ay' )
%     disp( Ay )
%     
%     disp( 'bx' )
%     disp( bx )
%     disp( 'by' )
%     disp( by )
%     
%     pause
    
    new_xg = Ax \ bx;
    new_yg = Ay \ by;
    
    xx_g( 2 : end - 1, 2 : end - 1 ) = reshape( new_xg, grid_sz_i, grid_sz_j );
    yy_g( 2 : end - 1, 2 : end - 1 ) = reshape( new_yg, grid_sz_i, grid_sz_j );
    
    
    clf
    %axis( [-1,1,-1,1] )
    axis auto
    hold on
    
    for s = 1 : 4
        plot( boundaries( s ).xx, boundaries( s ).yy, 'b' )
    end
    
    plotGrid( xx_g, yy_g )
    
    disp( 'iter: ' )
    disp( k )
    
    pause
end

% everything has an offset of 2, because of the boundary conditions

function g = c_alpha( xg, yg, i, j )
    g = 0.25 * ( ( xg( i + 1, j + 2 ) - xg( i + 1, j ) ) ^ 2 + ...
                 ( yg( i + 1, j + 2 ) - yg( i + 1, j ) ) ^ 2 );
end

function g = c_beta( xg, yg, i, j )
    g = 0.25 * ( ( xg( i + 2, j + 1 ) - xg( i, j + 1 ) ) * ...
                 ( xg( i + 1, j + 2 ) - xg( i + 1, j ) ) + ...
                 ( yg( i + 2, j + 1 ) - yg( i, j + 1 ) ) * ...
                 ( yg( i + 1, j + 2 ) - yg( i + 1, j ) ) );
end

function g = c_gamma( xg, yg, i, j )
    g = 0.25 * ( ( xg( i + 2, j + 1 ) - xg( i, j + 1 ) ) ^ 2 + ...
                 ( yg( i + 2, j + 1 ) - yg( i, j + 1 ) ) ^ 2 );
end

function indx = getIndxInMatrix( i, j, gSz_i, ~ )
    indx = i + ( j - 1 ) * gSz_i;
end

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

function valBoundary = getBoundaryConditionY( xg, yg, i, j, in, jn, gSz_i, gSz_j )

    di = in - i;
    dj = jn - j;

    if di == 1 && dj == 1

        valBoundary =  -0.5 * c_beta( xg, yg, i, j ) * yg( in + 1, jn + 1 );

    elseif di == 1 && dj == 0

        valBoundary =        c_alpha( xg, yg, i, j ) * yg( in + 1, jn + 1 );

    elseif di == 1 && dj == -1

        valBoundary =   0.5 * c_beta( xg, yg, i, j ) * yg( in + 1, jn + 1 );

    elseif di == 0 && dj == -1

        valBoundary =        c_gamma( xg, yg, i, j ) * yg( in + 1, jn + 1 );

    elseif di == -1 && dj == -1

        valBoundary = -0.5 * c_beta( xg, yg, i , j ) * yg( in + 1, jn + 1 );

    elseif di == -1 && dj == 0

        valBoundary =        c_alpha( xg, yg, i, j ) * yg( in + 1, jn + 1 );

    elseif di == -1 && dj == 1

        valBoundary =   0.5 * c_beta( xg, yg, i, j ) * yg( in + 1, jn + 1 );

    elseif di == 0 && dj == 1

        valBoundary =        c_gamma( xg, yg, i, j ) * yg( in + 1, jn + 1 );
        
    else 
        disp( 'shouldnt get here ' );
    end

end

%%%%%%%%%%%%%%%%%%%
%    
%    
%
%%%%%%%%%%%%%%%%%%%

function [A, b] = applyStencilX( A, b, xg, yg, i, j, gSz_i, gSz_j )

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
            
            b( s_row, 1 ) = -getBoundaryConditionX( xg, yg, i, j, i + 1, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i - 1, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i - 1, j, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i - 1, j + 1, gSz_i, gSz_j );
            
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
            
            b( s_row, 1 ) = -getBoundaryConditionX( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i    , j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i - 1, j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i - 1, j    , gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i - 1, j - 1, gSz_i, gSz_j );
            
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
            
            b( s_row, 1 ) = -getBoundaryConditionX( xg, yg, i, j, i - 1, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i - 1, j    , gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i - 1, j + 1, gSz_i, gSz_j );
            
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
            
            b( s_row, 1 ) = -getBoundaryConditionX( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i + 1, j    , gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i + 1, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i    , j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i - 1, j - 1, gSz_i, gSz_j );
            
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
            
            b( s_row, 1 ) = -getBoundaryConditionX( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i + 1, j    , gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i    , j + 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i - 1, j + 1, gSz_i, gSz_j );
            
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
            
            b( s_row, 1 ) = -getBoundaryConditionX( xg, yg, i, j, i + 1, j - 1, gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i + 1, j    , gSz_i, gSz_j ) - ...
                             getBoundaryConditionX( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j );            
            
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
        
        b( s_row, 1 ) = -getBoundaryConditionX( xg, yg, i, j, i + 1, j - 1, gSz_i, gSz_j ) - ...
                         getBoundaryConditionX( xg, yg, i, j, i    , j - 1, gSz_i, gSz_j ) - ...
                         getBoundaryConditionX( xg, yg, i, j, i - 1, j - 1, gSz_i, gSz_j );            
        
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
        
        b( s_row, 1 ) = -getBoundaryConditionX( xg, yg, i, j, i + 1, j + 1, gSz_i, gSz_j ) - ...
                         getBoundaryConditionX( xg, yg, i, j, i    , j + 1, gSz_i, gSz_j ) - ...
                         getBoundaryConditionX( xg, yg, i, j, i - 1, j + 1, gSz_i, gSz_j );            
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



function plotGrid( xg, yg )

    for i = 1 : size( xg, 1 ) 
       
        for j = 1 : size( xg, 2 )
           
            if i > 1
                plot( [ xg( i - 1, j ), xg( i, j ) ], ...
                      [ yg( i - 1, j ), yg( i, j ) ], 'k' );
            end
            if j > 1
                plot( [ xg( i, j - 1 ), xg( i, j ) ], ...
                      [ yg( i, j - 1 ), yg( i, j ) ], 'k' );
            end
            
            if i < size( xg, 1 )
                plot( [ xg( i + 1, j ), xg( i, j ) ], ...
                      [ yg( i + 1, j ), yg( i, j ) ], 'k' );
            end
            if j < size( xg, 2 )
                plot( [ xg( i, j + 1 ), xg( i, j ) ], ...
                      [ yg( i, j + 1 ), yg( i, j ) ], 'k' );
            end
            
        end
        
    end

end