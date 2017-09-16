clc
clear all
close all

xMin = -1;
xMax = 1;
yMin = -1;
yMax = 1;

boundaries(1).xx = []; boundaries(1).yy = [];
boundaries(2).xx = []; boundaries(2).yy = [];
boundaries(3).xx = []; boundaries(3).yy = [];
boundaries(4).xx = []; boundaries(4).yy = [];

option = input( 'geometry? (1:square-hole) (2:circle-hole) (3:circle) (4:square)' );

switch option
    
    case 1
        
        r = 0.5;
        
        N1 = 10;
        boundaries(1).xx = [ boundaries(1).xx, -1 * ones( 1, N1 ) ];
        boundaries(1).xx = [ boundaries(1).xx, linspace( -1, 0, N1 )];
        boundaries(1).yy = [ boundaries(1).yy, linspace( 0, -1, N1 ) ];
        boundaries(1).yy = [ boundaries(1).yy, -1 * ones( 1, N1 ) ];
        
        N2 = 20;
        boundaries(2).xx = [ boundaries(2).xx, zeros( 1, N2 ) ];
        boundaries(2).yy = [ boundaries(2).yy, linspace( -1, -r, N2 ) ];
        
        N3 = 20;
        theta = linspace( 1.5 * pi, pi, N3 );
        boundaries(3).xx = [ boundaries(3).xx, r * cos( theta ) ];
        boundaries(3).yy = [ boundaries(3).yy, r * sin( theta ) ];
    
        N4 = 20;
        boundaries(4).xx = [ boundaries(4).xx, linspace( -r, -1, N4 ) ];
        boundaries(4).yy = [ boundaries(4).yy, zeros( 1, N4 ) ];
        
    case 2
        
        r = 0.5;
        
        N1 = 20;
        theta = linspace( pi, 1.5 * pi, N1 );
        boundaries(1).xx = [ boundaries(1).xx, 2 * r * cos( theta ) ];
        boundaries(1).yy = [ boundaries(1).yy, 2 * r * sin( theta ) ];
        
        N2 = 20;
        boundaries(2).xx = [ boundaries(2).xx, zeros( 1, N2 ) ];
        boundaries(2).yy = [ boundaries(2).yy, linspace( -1, -r, N2 ) ];
        
        N3 = 20;
        theta = linspace( 1.5 * pi, pi, N3 );
        boundaries(3).xx = [ boundaries(3).xx, r * cos( theta ) ];
        boundaries(3).yy = [ boundaries(3).yy, r * sin( theta ) ];
    
        N4 = 20;
        boundaries(4).xx = [ boundaries(4).xx, linspace( -r, -1, N4 ) ];
        boundaries(4).yy = [ boundaries(4).yy, zeros( 1, N4 ) ];
        
    case 3
        
        r = 0.5;
        
        N = 20;
        theta = linspace( pi, 1.5 * pi, N );
        boundaries(1).xx = [ boundaries(1).xx, r * cos( theta ) ];
        boundaries(1).yy = [ boundaries(1).yy, r * sin( theta ) ];

        theta = linspace( 1.5 * pi, 2 * pi, N );
        boundaries(2).xx = [ boundaries(2).xx, r * cos( theta ) ];
        boundaries(2).yy = [ boundaries(2).yy, r * sin( theta ) ];
        
        theta = linspace( 0, 0.5 * pi, N );
        boundaries(3).xx = [ boundaries(3).xx, r * cos( theta ) ];
        boundaries(3).yy = [ boundaries(3).yy, r * sin( theta ) ];
        
        theta = linspace( 0.5 * pi, pi, N );
        boundaries(4).xx = [ boundaries(4).xx, r * cos( theta ) ];
        boundaries(4).yy = [ boundaries(4).yy, r * sin( theta ) ];
        
    case 4
        
        N = 20;
        boundaries(1).xx = [ boundaries(1).xx, zeros( 1, N ) ];
        boundaries(1).yy = [ boundaries(1).yy, linspace( 1, 0, N ) ];
        
        boundaries(2).xx = [ boundaries(2).xx, linspace( 0, 1, N ) ];
        boundaries(2).yy = [ boundaries(2).yy, zeros( 1, N ) ];
        
        boundaries(3).xx = [ boundaries(3).xx, ones( 1, N ) ];
        boundaries(3).yy = [ boundaries(3).yy, linspace( 0, 1, N ) ];
        
        boundaries(4).xx = [ boundaries(4).xx, linspace( 1, 0, N ) ];
        boundaries(4).yy = [ boundaries(4).yy, ones( 1, N ) ];
end

fg = figure(1);
axis( [ 2 * xMin, 2 * xMax, 2 * yMin, 2 * yMax ] );
hold on;

for q = 1 : 4
    plot( boundaries(q).xx, boundaries(q).yy );
end

fileHandle = fopen( "geometry_2_" + option + ".txt", 'w' );
for q = 1 : 4
    fprintf( fileHandle, '%f %f \n', [boundaries(q).xx; boundaries(q).yy] );
    fprintf( fileHandle, '---\n' );
end

fclose( fileHandle );


