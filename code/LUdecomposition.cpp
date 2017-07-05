
#include <cassert>
#include <iostream>
#include <armadillo>
#include <vector>

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

            LUcalc doolittle( const mat &A )
            {
                LUcalc _res;

                // Must implement these correspondences
                // U( i, j ) = ( A(i,j) - \sum_{k=1}^{i-1}L(i,k)U(k,j) )
                // L( i, j ) = \frac{ A(i,j) - \sum_{k=1}^{j-1}L(i,k)U(k,j) }{U(j,j)}

                assert( A.n_cols == A.n_rows );

                int N = A.n_cols;
                arma::mat _L = arma::eye<arma::mat>( A.n_rows, A.n_cols );
                arma::mat _U = A;

                for ( int i = 0; i < N; i++ )
                {
                    // At row i, calculate U( i, j ) for every j, as well as L( i, j )
                    for ( int j = 0; j < N; j++ )
                    {
                        if ( i <= j )
                        {
                            _U( i, j ) = A( i, j );
                            for ( int k = 0; k <= i - 1; k++ )
                            {
                                _U( i, j ) -= _L( i, k ) * _U( k, j );
                            }
                        }
                        else
                        {
                            _U( i, j ) = 0;
                        }

                        if ( i > j )
                        {
                            _L( i, j ) = A( i, j );
                            for ( int k = 0; k <= j - 1; k++ )
                            {
                                _L( i, j ) -= _L( i, k ) * _U( k, j );
                            }
                            _L( i, j ) = _L( i, j ) / _U( j, j );
                        }
                        else if ( i == j )
                        {
                            _L( i, j ) = 1;
                        }
                        else
                        {
                            _L( i, j ) = 0;
                        }
                    }
                }

                _res.L = _L;
                _res.U = _U;

                return _res;                
            }

            LUcalc crout( const mat &A )
            {
                LUcalc _res;

                // Must implement these correspondences
                // L(i,j) = ( A(i,j) - \sum_{k=1}^{j-1}L(i,k)U(k,j) )
                // U(i,j) = \frac{ A(i,j) - \sum_{k=1}^{i-1}L(i,k)U(k,j) }{L(i,i)}

                assert( A.n_cols == A.n_rows );

                int N = A.n_cols;
                arma::mat _U = arma::eye<arma::mat>( A.n_rows, A.n_cols );
                arma::mat _L = A;

                for ( int i = 0; i < N; i++ )
                {
                    // At row i, calculate L( i, j ) for every j, as well as U( i, j )
                    for ( int j = 0; j < N; j++ )
                    {
                        if ( i >= j )
                        {
                            _L( i, j ) = A( i, j );
                            for ( int k = 0; k <= j - 1; k++ )
                            {
                                _L( i, j ) -= _L( i, k ) * _U( k, j );
                            }
                        }
                        else
                        {
                            _L( i, j ) = 0;
                        }

                        if ( i < j )
                        {
                            _U( i, j ) = A( i, j );
                            for ( int k = 0; k <= i - 1; k++ )
                            {
                                _U( i, j ) -= _L( i, k ) * _U( k, j );
                            }
                            _U( i, j ) = _U( i, j ) / _L( i, i );
                        }
                        else if ( i == j )
                        {
                            _U( i, j ) = 1;
                        }
                        else
                        {
                            _U( i, j ) = 0;
                        }
                    }
                }

                _res.L = _L;
                _res.U = _U;

                return _res;
            }

            LUcalc gaussian( const mat &A )
            {
                assert( A.n_cols == A.n_rows );

                LUcalc _res;

                int N = A.n_cols;

                vector<double> _mk;
                for ( int q = 0; q < N; q++ )
                {
                    _mk.push_back( 0 );
                }

                // Set the initial state of the matrices to be calculated
                arma::mat _L = arma::eye<arma::mat>( A.n_rows, A.n_cols );
                arma::mat _U = A;
                arma::mat _Ak = A;
                arma::mat _Ek = arma::eye<arma::mat>( A.n_rows, A.n_cols );

                for ( int q = 0; q < N - 1; q++ )
                {
                    // initialize with 0s the mk vector of multipliers
                    for( int p = 0; p <= q; p++ )
                    {
                        _mk[p] = 0;
                    }

                    // initialize the other elements of mk with multipliers from Ak
                    for ( int p = q + 1; p < N; p++ )
                    {
                        // We are at column ( q + 1 ) ( actually q, for armadillo )
                        // extract the multipliers from col q of Ak
                        _mk[p] = -_Ak( p, q ) / _Ak( q, q );
                        //cout << "step" << endl;
                        //cout << _Ak( p, q ) << endl;
                        //cout << _Ak( q, q ) << endl;
                        //cout << _mk[p] << endl;

                        // populate L with this column of multipliers
                        // L = I - mk*ek^t, fill column k = q of L, from row p to N - 1
                        _L( p, q ) = -_mk[p];
                        
                        // Generate A(k+1) = Ek * Ak

                        // Eliminate the elements in column k accordingly
                        _Ak( p, q ) = 0;
                        for ( int c = q + 1; c < N; c++ )
                        {
                            _Ak( p, c ) += _mk[p] * _Ak( q, c );
                        }

                        //cout << _L << endl;
                    }
                }

                // U = A(N-1)
                _U = _Ak;

                _res.L = _L;
                _res.U = _U;

                return _res;
            }
            
        }
    }

}




int main()
{

    arma::mat A = { { 1,  5,  0,  0 },
                    { 2, 12,  5,  0 },
                    { 0,  4, 13,  5 },
                    { 0,  0,  6, 11 } };

    cout << "matrix A: " << endl;
    cout << A << endl;

    cout << "LU-Gaussian decomposition: --------------------------------" << endl;
    NLA::decompositions::lu::LUcalc _res = NLA::decompositions::lu::gaussian( A );
    cout << "L: " << endl;
    cout << _res.L << endl;
    cout << "U: " << endl;
    cout << _res.U << endl;

    cout << "LU-Crout decomposition: -----------------------------------" << endl;
    _res = NLA::decompositions::lu::crout( A );
    cout << "L: " << endl;
    cout << _res.L << endl;
    cout << "U: " << endl;
    cout << _res.U << endl;

    cout << "LU-Doolittle decomposition: -----------------------------------" << endl;
    _res = NLA::decompositions::lu::doolittle( A );
    cout << "L: " << endl;
    cout << _res.L << endl;
    cout << "U: " << endl;
    cout << _res.U << endl;


    return 0;
}