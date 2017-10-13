
clc
clear all
close all

xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;

fg = figure(1);
axis( [xMin, xMax, yMin, yMax ] );

global ax;
ax = axes;

global boundaryId;
boundaryId = 1;

global boundary1x; global boundary1y;
global boundary2x; global boundary2y;
global boundary3x; global boundary3y;
global boundary4x; global boundary4y;

boundary1x = []; boundary1y = [];
boundary2x = []; boundary2y = [];
boundary3x = []; boundary3y = [];
boundary4x = []; boundary4y = [];

plot( xMin, xMin, xMax, xMax, ...
      yMin, yMax, yMin, yMax );
hold on;

set( ax, 'buttondownfcn', @onClickCallback );
  
disp( 'fun' )
  
function onClickCallback( src, ev )
    global boundaryId;
    global boundary1x; global boundary1y;
    global boundary2x; global boundary2y;
    global boundary3x; global boundary3y;
    global boundary4x; global boundary4y;
    
    disp( 'foo' )
    %disp( ev )
    if ev.Button == 3
        % Change to next boundary
        disp( 'changed to next boundary' )
        boundaryId = boundaryId + 1;
        disp( boundaryId );
        
        if boundaryId == 5 
            set( src, 'buttondownfcn', @onDummyCallback );
            
            % Save current boundaries
            fileHandle = fopen( 'geometry_5.txt', 'w' );
%             fprintf( fileHandle, '%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n', ...
%                      boundary1x, boundary1y, ...
%                      boundary2x, boundary2y, ...
%                      boundary3x, boundary3y, ...
%                      boundary4x, boundary4y );
                 
            fprintf( fileHandle, '%f %f \n', [boundary1x; boundary1y] );
            fprintf( fileHandle, '---\n' );
            fprintf( fileHandle, '%f %f \n', [boundary2x; boundary2y] );
            fprintf( fileHandle, '---\n' );
            fprintf( fileHandle, '%f %f \n', [boundary3x; boundary3y] );
            fprintf( fileHandle, '---\n' );
            fprintf( fileHandle, '%f %f \n', [boundary4x; boundary4y] );
            fprintf( fileHandle, '---\n' );
            fclose( fileHandle );
        end
    else
        %%disp( ev.IntersectionPoint );
        x = ev.IntersectionPoint( 1, 1 );
        y = ev.IntersectionPoint( 1, 2 );
   
        %disp( x );
        %disp( y );
        
        color = 'bo';
        switch boundaryId 
            
            case 1
                color = 'bo';
                boundary1x( end + 1 ) = x;
                boundary1y( end + 1 ) = y;
            case 2
                color = 'ro';
                boundary2x( end + 1 ) = x;
                boundary2y( end + 1 ) = y;
            case 3
                color = 'go';
                boundary3x( end + 1 ) = x;
                boundary3y( end + 1 ) = y;
            case 4
                color = 'ko';
                boundary4x( end + 1 ) = x;
                boundary4y( end + 1 ) = y;
        end
        plot( src, x, y, color );
    end
end

function onDummyCallback( ~, ~ )
    disp( 'foo' );
end