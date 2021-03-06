
clear all
close all
clc

% define grid size in computational space
Ne = 30;
Nn = 30;

% create figures
fig1 = figure( 1 );
fig2 = figure( 2 );

% Load the geometry
boundaries = util_load_geometry( 'geometry_2_1.txt' );

% Use the algebraic generator
[xg, yg] = alg_transfinite_generator( boundaries, Ne, Nn );

% Plot the result of the algebraic generator
util_plot_grid( fig1, boundaries, xg, yg );

% for testing only
% [Ax, Ay, bx, by] = generateSysMatrix( xg, yg, size( xg, 1 ), size( xg, 2 ), size( xg, 1 ) - 2, size( xg, 2 ) - 2 );

% Define the number of iterations for the elliptic generator
nIters = 5;

% Use the elliptic generator

for q = 1 : nIters

	fprintf( 'solving iteration %d ...\n', q );

	[xg, yg] = elliptic_generator( xg, yg );

	figure( fig2 );
	axis auto
	hold on 

	clf;	
	util_plot_grid( fig2, boundaries, xg, yg );

	fprintf( 'solved iteration %d\n', q );

	pause

end