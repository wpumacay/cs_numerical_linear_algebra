
clc
clear all
close all

% Assume that the points are written to a file in a round like order ...
% like -> lastPoint_boundary_i = firstPoint_boundary_i+1

fileHandle = fopen( 'geometry_2_3.txt', 'r' );

% Read four boundaries
for q = 1 : 4

    boundaries( q ).xx = [];
    boundaries( q ).yy = [];
    
    while ( true )
        
        tline = fgetl( fileHandle );
        if strcmp( tline, '---' )
            % go to next boundary
            break
        end
        
        if ischar( tline )
            xy = str2num( tline );
            boundaries( q ).xx( end + 1 ) = xy( 1 );
            boundaries( q ).yy( end + 1 ) = xy( 2 );
        end
        
    end
    
end

% complete end points ( last_i = first_i+1 )
for q = 1 : 4
    b_size = size( boundaries( q ).xx );
    len_i = b_size( 2 );
    
    boundaries( q ).xx( 1, end + 1 ) = ...
        boundaries( rem( q , 4 ) + 1 ).xx( 1, 1 );
    
    boundaries( q ).yy( 1, end + 1 ) = ...
        boundaries( rem( q , 4 ) + 1 ).yy( 1, 1 );
    
    boundaries( q ).size = len_i + 1;
end

% correct the direction of the boundaries 
% left must go from bottom to top
% top must go from left to right
boundaries( 1 ).xx = flip( boundaries( 1 ).xx );
boundaries( 1 ).yy = flip( boundaries( 1 ).yy );

boundaries( 4 ).xx = flip( boundaries( 4 ).xx );
boundaries( 4 ).yy = flip( boundaries( 4 ).yy );

disp( boundaries );


% define each boundary using piecewise segments
for q = 1 : 4
    boundaries( q ).x_geo = [];
    boundaries( q ).y_geo = [];
    
    % calculate udir and len
    for i = 1 : ( boundaries( q ).size - 1 )
        x_p0 = boundaries( q ).xx( i );
        x_p1 = boundaries( q ).xx( i + 1 );
        y_p0 = boundaries( q ).yy( i );
        y_p1 = boundaries( q ).yy( i + 1 );
        
        dx = x_p1 - x_p0; dy = y_p1 - y_p0;
        len = sqrt( dx^2 + dy^2 );
        
        ux = dx / len; uy = dy / len;
        l_geo_x = [x_p0, ux, len];
        l_geo_y = [y_p0, uy, len];
        
        boundaries( q ).x_geo( end + 1, : ) = l_geo_x;
        boundaries( q ).y_geo( end + 1, : ) = l_geo_y;
    end
end

% define grid dimensions ( xi=e, eta=n )
Ne = 20;
Nn = 20;

de = 1 / Ne;
dn = 1 / Nn;

ee = linspace( 0, 1, Ne + 1 );
nn = linspace( 0, 1, Nn + 1 );

xx_g = zeros( Nn - 1, Ne - 1 );
yy_g = zeros( Nn - 1, Ne - 1 );

for i = 2 : Nn
    
   for j = 2 : Ne 
      
       e = ee( j );
       n = nn( i );
       
       indxs_e = zeros( 1, 4 );
       indxs_n = zeros( 1, 4 );
       
       % calculate the index of which segment to use for each boundary
       for s = 1 : 4
           n_segments = boundaries( s ).size - 1;
           b_d = 1 / n_segments;
           
           indxs_e( s ) = floor( e / b_d ) + 1;
           indxs_n( s ) = floor( n / b_d ) + 1;
       end
      
       xl = 0; yl = 0;
       xr = 0; yr = 0;
       xb = 0; yb = 0;
       xt = 0; yt = 0;
       
       xb_0 = boundaries( 2 ).xx( 1 ); yb_0 = boundaries( 2 ).yy( 1 );
       xb_1 = boundaries( 2 ).xx( end ); yb_1 = boundaries( 2 ).yy( end );
       
       xt_0 = boundaries( 4 ).xx( 1 ); yt_0 = boundaries( 4 ).yy( 1 );
       xt_1 = boundaries( 4 ).xx( end ); yt_1 = boundaries( 4 ).yy( end );
       
       % calculate xl, xr, xb and xt using the appropiate boundaries
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
       
       x = ( 1 - e ) * xl + ( e ) * xr + ( 1 - n ) * xb + ( n ) * xt - ...
           ( 1 - n ) * ( 1 - e ) * xb_0 - ( 1 - e ) * n * xt_0 - ...
           ( 1 - n ) * e * xb_1 - n * e * xt_1;

       y = ( 1 - e ) * yl + ( e ) * yr + ( 1 - n ) * yb + ( n ) * yt - ...
           ( 1 - n ) * ( 1 - e ) * yb_0 - ( 1 - e ) * n * yt_0 - ...
           ( 1 - n ) * e * yb_1 - n * e * yt_1;
       
       xx_g( i - 1, j - 1 ) = x; yy_g( i - 1, j - 1 ) = y;
       
   end
    
end

figure(1)
axis( [-2,2,-2,2] )
hold on

for s = 1 : 4
    plot( boundaries( s ).xx, boundaries( s ).yy, 'b' )
end

for i = 1 : size( xx_g, 1 )
   
    for j = 1 : size( xx_g, 2 )
       
        plot( xx_g( i, j ) , yy_g( i, j ), 'k.' )
        
    end
    
end

plotGrid( xx_g, yy_g ) 

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


