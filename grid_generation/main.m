
clear all
close all
clc

% define grid size in computational space
Ne = 20;
Nn = 20;

% create figures
fig1 = figure( 1 );
fig2 = figure( 2 );

% Load the geometry
boundaries = util_load_geometry( 'geometry_2_1.txt' );

% Use the algebraic generator
[xg, yg] = alg_transfinite_generator( boundaries, Ne, Nn );

% Plot the result of the algebraic generator
util_plot_grid( fig1, boundaries, xg, yg );

% Define the number of iterations for the elliptic generator
nIters = 10;

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

	pause;

end