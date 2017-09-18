
% @brief: generates a grid using the algebraic method based on transfinite interpolation
% @param {array[struct]} boundaries : boundaries of the grid
% @param {int} Ne : size of the grid in curvilinear coordinates ( e (xi) coordinate )
% @param {int} Nn : size of the grid in curvilinear coordinates ( n (eta) coordinate )
% @return {[matrix, matrix]} [xx_g, yy_g] : returns the x and y matrices representing the generated grid
function [xx_g, yy_g] = alg_transfinite_generator( boundaries, Ne, Nn )

	% create grid in curvilinear coordinates %%%%%%%
	ee = linspace( 0, 1, Ne + 1 );
	nn = linspace( 0, 1, Nn + 1 );
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% create the grid in physical coordinates x and y %%%%%%%%
	xx_g = zeros( Nn + 1, Ne + 1 );
	yy_g = zeros( Nn + 1, Ne + 1 );
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Transfinite interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 'e' (xi)  relates to the 'i' direction
	% 'n' (eta) relates to the 'j' direction
	for i = 1 : ( Ne + 1 )
	    
	   for j = 1 : ( Nn + 1 )
	       
	       % Get e,n coordinates in computational space %%%%%
	       e = ee( i );
	       n = nn( j );
	       % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	       
	       indxs_e = zeros( 1, 4 );
	       indxs_n = zeros( 1, 4 );
	       
	       % calculate the index of which segment to use respect to each boundary %%%%%%
	       for s = 1 : 4
	           n_segments = boundaries( s ).size - 1;
	           b_d = 1 / n_segments;
	           
	           indxs_e( s ) = floor( e / b_d ) + 1;
	           if indxs_e( s ) > n_segments
	               indxs_e( s ) = n_segments;
	           end
	           
	           indxs_n( s ) = floor( n / b_d ) + 1;
	           if indxs_n( s ) > n_segments
	               indxs_n( s ) = n_segments;
	           end
	       end
	       % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	      
	       % Initialize variables for transfinite interpolatino %%%%%%%%%%%%%%%%%%%%%%%%
	       xl = 0; yl = 0;
	       xr = 0; yr = 0;
	       xb = 0; yb = 0;
	       xt = 0; yt = 0;
	       
	       xb_0 = boundaries( 2 ).xx( 1 ); yb_0 = boundaries( 2 ).yy( 1 );
	       xb_1 = boundaries( 2 ).xx( end ); yb_1 = boundaries( 2 ).yy( end );
	       
	       xt_0 = boundaries( 4 ).xx( 1 ); yt_0 = boundaries( 4 ).yy( 1 );
	       xt_1 = boundaries( 4 ).xx( end ); yt_1 = boundaries( 4 ).yy( end );
	       % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	       
	       % calculate xl, xr, xb and xt using the appropiate boundaries %%%%%%%%%%%%%%%
	       for s = 1 : 4
	           
	           n_segments = boundaries( s ).size - 1;
	           b_d = 1 / n_segments;
	           
	           switch s
	               
	               case 1
	                   seg_params_x = boundaries( s ).x_geo( indxs_n( s ), : );
	                   seg_params_y = boundaries( s ).y_geo( indxs_n( s ), : );
	                   
	                   x_p0 = seg_params_x( 1 ); y_p0 = seg_params_y( 1 );
	                   udx = seg_params_x( 2 ); udy = seg_params_y( 2 );
	                   lenx = seg_params_x( 3 ); leny = seg_params_y( 3 );
	                   
	                   xl = x_p0 + ( n - ( indxs_n( s ) - 1 ) * b_d ) * n_segments * udx * lenx;
	                   yl = y_p0 + ( n - ( indxs_n( s ) - 1 ) * b_d ) * n_segments * udy * leny;
	                   
	               case 2
	                   seg_params_x = boundaries( s ).x_geo( indxs_e( s ), : );
	                   seg_params_y = boundaries( s ).y_geo( indxs_e( s ), : );
	                   
	                   x_p0 = seg_params_x( 1 ); y_p0 = seg_params_y( 1 );
	                   udx = seg_params_x( 2 ); udy = seg_params_y( 2 );
	                   lenx = seg_params_x( 3 ); leny = seg_params_y( 3 );
	                   
	                   xb = x_p0 + ( e - ( indxs_e( s ) - 1 ) * b_d ) * n_segments * udx * lenx;
	                   yb = y_p0 + ( e - ( indxs_e( s ) - 1 ) * b_d ) * n_segments * udy * leny;
	                   
	               case 3
	                   seg_params_x = boundaries( s ).x_geo( indxs_n( s ), : );
	                   seg_params_y = boundaries( s ).y_geo( indxs_n( s ), : );
	                   
	                   x_p0 = seg_params_x( 1 ); y_p0 = seg_params_y( 1 );
	                   udx = seg_params_x( 2 ); udy = seg_params_y( 2 );
	                   lenx = seg_params_x( 3 ); leny = seg_params_y( 3 );
	                   
	                   xr = x_p0 + ( n - ( indxs_n( s ) - 1 ) * b_d ) * n_segments * udx * lenx;
	                   yr = y_p0 + ( n - ( indxs_n( s ) - 1 ) * b_d ) * n_segments * udy * leny;
	                   
	               case 4
	                   seg_params_x = boundaries( s ).x_geo( indxs_e( s ), : );
	                   seg_params_y = boundaries( s ).y_geo( indxs_e( s ), : );
	                   
	                   x_p0 = seg_params_x( 1 ); y_p0 = seg_params_y( 1 );
	                   udx = seg_params_x( 2 ); udy = seg_params_y( 2 );
	                   lenx = seg_params_x( 3 ); leny = seg_params_y( 3 );
	                   
	                   xt = x_p0 + ( e - ( indxs_e( s ) - 1 ) * b_d ) * n_segments * udx * lenx;
	                   yt = y_p0 + ( e - ( indxs_e( s ) - 1 ) * b_d ) * n_segments * udy * leny;
	           end
	           
	       end
	       % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	       % Apply transfinite interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	       x = ( 1 - e ) * xl + ( e ) * xr + ( 1 - n ) * xb + ( n ) * xt - ...
	           ( 1 - n ) * ( 1 - e ) * xb_0 - ( 1 - e ) * n * xt_0 - ...
	           ( 1 - n ) * e * xb_1 - n * e * xt_1;

	       y = ( 1 - e ) * yl + ( e ) * yr + ( 1 - n ) * yb + ( n ) * yt - ...
	           ( 1 - n ) * ( 1 - e ) * yb_0 - ( 1 - e ) * n * yt_0 - ...
	           ( 1 - n ) * e * yb_1 - n * e * yt_1;
	       
	       xx_g( i, j ) = x; yy_g( i, j ) = y;
	       % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	       
	   end
	    
	end
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end