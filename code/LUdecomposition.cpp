
#include <cassert>
#include <iostream>
#include <armadillo>
#include <vector>
#include <cmath>

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
                mat P;
                mat Q;
            };

            namespace pivoting
            {
                enum _pivoting
                {
                    NO_PIVOTING,
                    PARTIAL_PIVOTING,
                    FULL_PIVOTING
                };
            }

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

            LUcalc gaussian( const mat &A, pivoting::_pivoting pPivoting = pivoting::NO_PIVOTING )
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
                arma::mat _L = arma::zeros<arma::mat>( A.n_rows, A.n_cols );
                arma::mat _U = A;
                arma::mat _P = arma::eye<arma::mat>( A.n_rows, A.n_cols );
                arma::mat _Q = arma::eye<arma::mat>( A.n_rows, A.n_cols );

                arma::mat _Ak = A;
                arma::mat _Ek = arma::eye<arma::mat>( A.n_rows, A.n_cols );

                for ( int q = 0; q < N - 1; q++ )
                {
                    // initialize with 0s the mk vector of multipliers
                    for( int p = 0; p <= q; p++ )
                    {
                        _mk[p] = 0;
                    }

                    int _row_pivot_indx = q;// First, element q,q
                    if ( pPivoting == pivoting::PARTIAL_PIVOTING )
                    {
                        // Find the index of the biggest pivot
                        for ( int p = q + 1; p < N; p++ )
                        {
                            if ( abs( _Ak( _row_pivot_indx, q ) ) < abs( _Ak( p, q ) ) )
                            {
                                _row_pivot_indx = p;
                            }
                        }

                        if ( _Ak( _row_pivot_indx, q ) == 0 )
                        {
                            cout << "pivot is 0, danger!" << endl;
                        }
                        else
                        {
                            cout << "changed: " << _row_pivot_indx << " - " << q << endl;
                        }

                        // exchange rows in Ak and P
                        for ( int s = 0; s < N; s++ )
                        {
                            if ( s >= q )
                            {
                                double _tmp = _Ak( _row_pivot_indx, s );
                                _Ak( _row_pivot_indx, s ) = _Ak( q, s );
                                _Ak( q, s ) = _tmp;
                            }

                            double _tmp = _P( _row_pivot_indx, s );
                            _P( _row_pivot_indx, s ) = _P( q, s );
                            _P( q, s ) = _tmp;
                        }

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

                _res.L = _L + arma::eye<arma::mat>( N, N );
                _res.U = _U;
                _res.P = _P;
                _res.Q = _Q;

                return _res;
            }

            arma::mat cholesky( const mat &A )
            {
                assert( A.n_cols == A.n_rows );

                int N = A.n_cols;

                arma::mat _H = arma::zeros<mat>( N, N );

                for ( int p = 0; p < N; p++ )
                {
                    double _sum = 0;
                    for ( int k = 0; k <= p - 1; k++ )
                    {
                        _sum += ( _H( p, k ) * _H( p, k ) );
                    }
                    _H( p, p ) = sqrt( A( p, p ) - _sum );

                    for ( int q = p + 1; q < N; q++ )
                    {
                        double _sum = 0;
                        for ( int k = 0; k <= p - 1; k++ )
                        {
                            _sum += ( _H( q, k ) * _H( p, k ) );
                        }
                        _H( q, p ) = ( A( q, p ) - _sum ) / _H( p, p );
                    }
                }

                return _H;
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

    cout << "LU-Gaussian decomposition: ------------------------------------" << endl;
    NLA::decompositions::lu::LUcalc _res = NLA::decompositions::lu::gaussian( A );
    cout << "L: " << endl;
    cout << _res.L << endl;
    cout << "U: " << endl;
    cout << _res.U << endl;

    cout << "LU-Gaussian decomposition / partial pivoting: -----------------" << endl;
    //arma::mat S = { { 0,  1,  1 },
    //                { 2,  2,  3 },
    //                { 1,  2,  1 } };
    arma::mat S = { { 2,  1,  1,  0 },
                    { 4,  3,  3,  1 },
                    { 8,  7,  9,  5 },
                    { 6,  7,  9,  8 } };
    _res = NLA::decompositions::lu::gaussian( S, NLA::decompositions::lu::pivoting::PARTIAL_PIVOTING );
    cout << "matrix: " << endl;
    cout << S << endl;
    cout << "L: " << endl;
    cout << _res.L << endl;
    cout << "U: " << endl;
    cout << _res.U << endl;
    cout << "P: " << endl;
    cout << _res.P << endl;

    cout << "check: " << endl;
    cout << "PS: " << endl;
    cout << _res.P * S << endl;
    cout << "LU: " << endl;
    cout << _res.L * _res.U << endl;

    cout << "LU-Crout decomposition: ---------------------------------------" << endl;
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

    cout << "Cholesky decomposition: ---------------------------------------" << endl;
    arma::mat B = arma::randu<arma::mat>( 4, 4 );
    arma::mat C = B.t() * B;
    arma::mat _h = NLA::decompositions::lu::cholesky( C );
    cout << "H: " << endl;
    cout << _h << endl;
    cout << "H-arma: " << endl;
    cout << arma::chol( C, "lower" ) << endl;


    return 0;
}