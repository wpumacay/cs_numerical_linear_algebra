
% @brief: plots the geometry and a generated grid on top
% @param {figure} fig : figure handle
% @param {array[struct]} boundaries : array of structs that define the boundaries
% @param {matrix} xg : x-grid
% @param {matrix} yg : y-grid
function util_plot_grid( fig, boundaries, xg, yg )

	figure( fig ) % make current
	axis auto
	hold on

	% plot the given grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % plot the boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for s = 1 : 4
	    plot( boundaries( s ).xx, boundaries( s ).yy, 'b' )
	end
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



