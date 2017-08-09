
#include <iostream>
#include <armadillo>
#include <cmath>

#define PI 3.14159265358979

using namespace std;
using namespace arma;


namespace plate2Dsolver
{


    class Plate2D
    {
        // Discretization part
        double dx;
        double dy;
        int nx;
        int ny;

        double lx;
        double ly;

        arma::mat sysMat;
        arma::mat sysVec;
        arma::mat sysSol;
        int dim;

    public :

        Plate2D( double lenX, double lenY, double pDx, double pDy )
        {
            lx = lenX;
            ly = lenY;
            dx = pDx;
            dy = pDy;

            nx = ( ( int ) lx / dx ) + 1;// discrete points in x from 0 to nx - 1
            ny = ( ( int ) ly / dy ) + 1;// discrete points in y from 0 to ny - 1

            dim = ( ny - 2 ) * ( nx - 2 );
            sysMat = arma::zeros<arma::mat>( dim, dim );
            sysVec = arma::zeros<arma::mat>( dim, 1 );
            sysSol = arma::zeros<arma::mat>( dim, 1 );
        }       

        int getIndxInSysMat( int i, int j )
        {
            return  ( i - 1 ) * ( nx - 2 ) + ( j - 1 );
        }

        // Modify this accordingly. In this case assumed a square of len 1 as plate ...
        // and initial conditions as described in a problem in a book
        double getBoundaryCondition( int i, int j )
        {
            if ( i == 0 )
            {
                return 1.0 - j * dx;
            }
            else if ( i == ny - 1 )
            {
                return 1.0;
            }
            else if ( j == 0 )
            {
                return 1.0;
            }
            else if ( j == nx - 1 )
            {
                return i * dy;
            }
            cout << "shouldnt get boundary condition for : " << i << " - " << j << endl;
            return 0.0;
        }

        double f( int i, int j )
        {
            double x = j * dx;
            double y = i * dy;
            double _f = -PI * PI * sin( PI * x ) * sin( PI * y );
            return _f * dx * dx;
        }

        void applyStencil( int i, int j )
        {
            // Generate the corresponding row

            // Apply reduced stencil to generate the row of the matrix
            // because the i or j if extreme are boundary conditions

            int s_row = getIndxInSysMat( i, j );
            if ( i == 1 )
            {
                if ( j == 1 )
                {
                    // bottom left corner
                    int _center_col = s_row;
                    int _right_col  = getIndxInSysMat( i, j + 1 );
                    int _up_col     = getIndxInSysMat( i + 1, j );

                    sysMat( s_row, _center_col ) = -4;
                    sysMat( s_row, _right_col )  = 1;
                    sysMat( s_row, _up_col )     = 1;

                    sysVec( s_row, 0 ) = f( i, j ) - 
                                         getBoundaryCondition( i - 1, j ) -
                                         getBoundaryCondition( i, j - 1 );

                }
                else if ( j == nx - 2 )
                {
                    // bottom left corner
                    int _center_col = s_row;
                    int _left_col   = getIndxInSysMat( i, j - 1 );
                    int _up_col     = getIndxInSysMat( i + 1, j );

                    sysMat( s_row, _center_col ) = -4;
                    sysMat( s_row, _left_col )   = 1;
                    sysMat( s_row, _up_col )     = 1;

                    sysVec( s_row, 0 ) = f( i, j ) - 
                                         getBoundaryCondition( i - 1, j ) -
                                         getBoundaryCondition( i, j + 1 );
                }
                else
                {
                    // bottom side
                    int _center_col = s_row;
                    int _left_col   = getIndxInSysMat( i, j - 1 );
                    int _right_col  = getIndxInSysMat( i, j + 1 );
                    int _up_col     = getIndxInSysMat( i + 1, j );

                    sysMat( s_row, _center_col ) = -4;
                    sysMat( s_row, _left_col )   = 1;
                    sysMat( s_row, _right_col )  = 1;
                    sysMat( s_row, _up_col )     = 1;

                    sysVec( s_row, 0 ) = f( i, j ) - 
                                         getBoundaryCondition( i - 1, j );
                   
                }
            }
            else if ( i == ny - 2 )
            {
                if ( j == 1 )
                {
                    // top left corner
                    int _center_col = s_row;
                    int _right_col  = getIndxInSysMat( i, j + 1 );
                    int _down_col   = getIndxInSysMat( i - 1, j );

                    sysMat( s_row, _center_col ) = -4;
                    sysMat( s_row, _right_col )  = 1;
                    sysMat( s_row, _down_col )   = 1;

                    sysVec( s_row, 0 ) = f( i, j ) - 
                                         getBoundaryCondition( i + 1, j ) -
                                         getBoundaryCondition( i, j - 1 );

                }
                else if ( j == nx - 2 )
                {
                    // top right corner
                    int _center_col = s_row;
                    int _left_col   = getIndxInSysMat( i, j - 1 );
                    int _down_col   = getIndxInSysMat( i - 1, j );

                    sysMat( s_row, _center_col ) = -4;
                    sysMat( s_row, _left_col )   = 1;
                    sysMat( s_row, _down_col )   = 1;

                    sysVec( s_row, 0 ) = f( i, j ) - 
                                         getBoundaryCondition( i + 1, j ) -
                                         getBoundaryCondition( i, j + 1 );
                }
                else
                {
                    // up side
                    int _center_col = s_row;
                    int _left_col   = getIndxInSysMat( i, j - 1 );
                    int _right_col  = getIndxInSysMat( i, j + 1 );
                    int _down_col   = getIndxInSysMat( i - 1, j );

                    sysMat( s_row, _center_col ) = -4;
                    sysMat( s_row, _left_col )   = 1;
                    sysMat( s_row, _right_col )  = 1;
                    sysMat( s_row, _down_col )   = 1;

                    sysVec( s_row, 0 ) = f( i, j ) - 
                                         getBoundaryCondition( i + 1, j );
                }
            }
            else if ( j == 1 )
            {
                // left side
                int _center_col = s_row;
                int _up_col     = getIndxInSysMat( i + 1, j );
                int _right_col  = getIndxInSysMat( i, j + 1 );
                int _down_col   = getIndxInSysMat( i - 1, j );

                sysMat( s_row, _center_col ) = -4;
                sysMat( s_row, _up_col )     = 1;
                sysMat( s_row, _right_col )  = 1;
                sysMat( s_row, _down_col )   = 1;

                sysVec( s_row, 0 ) = f( i, j ) - 
                                     getBoundaryCondition( i, j - 1 );
           
            } 
            else if ( j == nx - 2 )
            {
                // right side
                int _center_col = s_row;
                int _left_col   = getIndxInSysMat( i, j - 1 );
                int _up_col     = getIndxInSysMat( i + 1, j );
                int _down_col   = getIndxInSysMat( i - 1, j );

                sysMat( s_row, _center_col ) = -4;
                sysMat( s_row, _left_col )   = 1;
                sysMat( s_row, _up_col )     = 1;
                sysMat( s_row, _down_col )   = 1;

                sysVec( s_row, 0 ) = f( i, j ) - 
                                     getBoundaryCondition( i, j + 1 );

            }
            else
            {
                // Apply the full stencil

                int _center_col = s_row;
                int _left_col   = getIndxInSysMat( i, j - 1 );
                int _right_col  = getIndxInSysMat( i, j + 1 );
                int _up_col     = getIndxInSysMat( i + 1, j );
                int _down_col   = getIndxInSysMat( i - 1, j );

                sysMat( s_row, _center_col ) = -4;
                sysMat( s_row, _left_col ) = 1;
                sysMat( s_row, _right_col ) = 1;
                sysMat( s_row, _up_col ) = 1;
                sysMat( s_row, _down_col ) = 1;

                sysVec( s_row, 0 ) = f( i, j );
            }
        }

        void generateGrid()
        {
            cout << "generating system matrix" << endl;
            // Generate the system matrix and system vector
            for ( int i = 1; i < ny - 1; i++ )
            {
                for ( int j = 1; j < nx - 1; j++ )
                {
                    // Apply stencil accordingly
                    applyStencil( i, j );
                }
            }

            

        }

        void solve()
        {
            generateGrid();
            
            sysSol = arma::solve( sysMat, sysVec );

            cout << "nx,ny: " << nx << " " << ny << endl;
            cout << "dim: " << dim << endl;

            cout << "sysMat *****" << endl;
            cout << sysMat << endl;
            cout << "************" << endl;

            cout << "sysVec *****" << endl;
            cout << sysVec << endl;
            cout << "************" << endl;

            cout << "sysSol ***" << endl;
            cout << sysSol << endl;
            cout << "************" << endl;

            cout << "actual sol**" << endl;
            for ( int i = 1; i < ny - 1; i++ )
            {
                for ( int j = 1; j < nx - 1; j++ )
                {
                    double x = j * dx;
                    double y = i * dy;
                    
                    double sol = 1 - x + x * y + 0.5 * sin( PI * x ) * sin( PI * y );
                    cout << sol << endl;
                }
            }
            cout << "************" << endl;
        }


    };









}
