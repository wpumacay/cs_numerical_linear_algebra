
disp( 'running alg. generator for 2x2 grid' );

alg_geometry_loader;

disp( 'finished running alg.generator for 2x2 grid' );

% calculate coefficients

alpha_11 = 0.25 * ( ( xx_g( 2, 3 ) - xx_g( 2, 1 ) )^2 + ...
                    ( yy_g( 2, 3 ) - yy_g( 2, 1 ) )^2 );
                
alpha_21 = 0.25 * ( ( xx_g( 3, 3 ) - xx_g( 3, 1 ) )^2 + ...
                    ( yy_g( 3, 3 ) - yy_g( 3, 1 ) )^2 );
                
alpha_12 = 0.25 * ( ( xx_g( 2, 4 ) - xx_g( 2, 2 ) )^2 + ...
                    ( yy_g( 2, 4 ) - yy_g( 2, 2 ) )^2 );
                
alpha_22 = 0.25 * ( ( xx_g( 3, 4 ) - xx_g( 3, 2 ) )^2 + ...
                    ( yy_g( 3, 4 ) - yy_g( 3, 2 ) )^2 );
                

gamma_11 = 0.25 * ( ( xx_g( 3, 2 ) - xx_g( 1, 2 ) )^2 + ...
                    ( yy_g( 3, 2 ) - yy_g( 1, 2 ) )^2 );
                
gamma_21 = 0.25 * ( ( xx_g( 4, 2 ) - xx_g( 2, 2 ) )^2 + ...
                    ( yy_g( 4, 2 ) - yy_g( 2, 2 ) )^2 );
                
gamma_12 = 0.25 * ( ( xx_g( 3, 3 ) - xx_g( 1, 3 ) )^2 + ...
                    ( yy_g( 3, 3 ) - yy_g( 1, 3 ) )^2 );
                
gamma_22 = 0.25 * ( ( xx_g( 4, 3 ) - xx_g( 2, 3 ) )^2 + ...
                    ( yy_g( 4, 3 ) - yy_g( 2, 3 ) )^2 );
                

beta_11 = 0.25 * ( ( xx_g( 3, 2 ) - xx_g( 1, 2 ) ) * ( xx_g( 2, 3 ) - xx_g( 2, 1 ) ) + ...
                   ( yy_g( 3, 2 ) - yy_g( 1, 2 ) ) * ( yy_g( 2, 3 ) - yy_g( 2, 1 ) ) );

beta_21 = 0.25 * ( ( xx_g( 4, 2 ) - xx_g( 2, 2 ) ) * ( xx_g( 3, 3 ) - xx_g( 3, 1 ) ) + ...
                   ( yy_g( 4, 2 ) - yy_g( 2, 2 ) ) * ( yy_g( 3, 3 ) - yy_g( 3, 1 ) ) );
               
beta_12 = 0.25 * ( ( xx_g( 3, 3 ) - xx_g( 1, 3 ) ) * ( xx_g( 2, 4 ) - xx_g( 2, 2 ) ) + ...
                   ( yy_g( 3, 3 ) - yy_g( 1, 3 ) ) * ( yy_g( 2, 4 ) - yy_g( 2, 2 ) ) );
               
beta_22 = 0.25 * ( ( xx_g( 4, 3 ) - xx_g( 2, 3 ) ) * ( xx_g( 3, 4 ) - xx_g( 3, 2 ) ) + ...
                   ( yy_g( 4, 3 ) - yy_g( 2, 3 ) ) * ( yy_g( 3, 4 ) - yy_g( 3, 2 ) ) );


Ax = [ ( -2 * alpha_11 - 2 * gamma_11 ) ,              alpha_11             ,               gamma_11            ,        -0.5 * beta_11             ;
                  alpha_21              , ( -2 * alpha_21 - 2 * gamma_21 )  ,           0.5 * beta_21           ,            gamma_21               ;
                  gamma_12              ,            0.5 * beta_12          , ( -2 * alpha_12 - 2 * gamma_12 )  ,            alpha_12               ;
               -0.5 * beta_22           ,              gamma_22             ,               alpha_22            , ( -2 * alpha_22 - 2 * gamma_22 )  ];
           
Ay = Ax;

bx = [ -alpha_11 * xx_g( 1, 2 ) - 0.5 * beta_11 * xx_g( 3, 1 ) - 0.5 * beta_11 * xx_g( 1, 3 ) + 0.5 * beta_11 * xx_g( 1, 1 ) - gamma_11 * xx_g( 2, 1 ) ;
       0.5 * beta_21 * xx_g( 4, 3 ) - alpha_21 * xx_g( 4, 2 ) - 0.5 * beta_21 * xx_g( 4, 1 ) - gamma_21 * xx_g( 3, 1 ) + 0.5 * beta_21 * xx_g( 2, 1 )  ;
       -alpha_12 * xx_g( 1, 3 ) - gamma_12 * xx_g( 2, 4 ) + 0.5 * beta_12 * xx_g( 3, 4 ) + 0.5 * beta_12 * xx_g( 1, 2 ) - 0.5 * beta_12 * xx_g( 1, 4 ) ;
       -alpha_22 * xx_g( 4, 3 ) - gamma_22 * xx_g( 3, 4 ) - 0.5 * beta_22 * xx_g( 2, 4 ) - 0.5 * beta_22 * xx_g( 4, 2 ) + 0.5 * beta_22 * xx_g( 4, 4 ) ];
                
by = [ -alpha_11 * yy_g( 1, 2 ) - 0.5 * beta_11 * yy_g( 3, 1 ) - 0.5 * beta_11 * yy_g( 1, 3 ) + 0.5 * beta_11 * yy_g( 1, 1 ) - gamma_11 * yy_g( 2, 1 ) ;
       0.5 * beta_21 * yy_g( 4, 3 ) - alpha_21 * yy_g( 4, 2 ) - 0.5 * beta_21 * yy_g( 4, 1 ) - gamma_21 * yy_g( 3, 1 ) + 0.5 * beta_21 * yy_g( 2, 1 )  ;
       -alpha_12 * yy_g( 1, 3 ) - gamma_12 * yy_g( 2, 4 ) + 0.5 * beta_12 * yy_g( 3, 4 ) + 0.5 * beta_12 * yy_g( 1, 2 ) - 0.5 * beta_12 * yy_g( 1, 4 ) ;
       -alpha_22 * yy_g( 4, 3 ) - gamma_22 * yy_g( 3, 4 ) - 0.5 * beta_22 * yy_g( 2, 4 ) - 0.5 * beta_22 * yy_g( 4, 2 ) + 0.5 * beta_22 * yy_g( 4, 4 ) ];
   

                
                
                
                
                
                
                
                
                
                
                
                
                
                