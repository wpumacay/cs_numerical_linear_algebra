
#include <cassert>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

namespace NLA
{


    namespace decompositions
    {

        namespace lu
        {
            
            struct LUcalc
            {
                mat L;
                mat U;
            };

            void crout( const mat A, mat &L, mat &U )
            {
                
            }

            LUcalc doolitle( const mat A )
            {
                LUcalc _res;
                _res.L = arma::eye<mat>( A.size() );
                _res.U = A;

                for ( int q = 0; q < A.n_cols; q++ )
                {

                }
            }
            
        }





    }


}




int main()
{





    return 0;
}