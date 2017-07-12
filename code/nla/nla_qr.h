

#pragma once

#include "nla.h"

using namespace std;
using namespace arma;

#define Sign( x ) ( x > 0 ? 1 : ( x == 0 ? 1 : -1 ) )

namespace NLA
{
    namespace decompositions
    {
        namespace qr
        {

            struct QRcalc
            {
                arma::mat Q;
                arma::mat R;
            };

            struct ZeroHcalc
            {
                arma::mat u;
                double sigma;
            };


            namespace householder
            {

                arma::mat mult_Hv( const arma::mat &u, const arma::mat &v )
                {
                    assert( u.n_cols == 1 && v.n_cols == 1 );

                    arma::mat _res;

                    // Hx = ( I - 2uu^{T}/u^{T}u )x = x - ( 2/u^{T}u )u(u^{T}x) = x - beta u u^{T}x
                    arma::vec _u = arma::vectorise( u );
                    arma::vec _v = arma::vectorise( v );

                    double _unorm = arma::norm( _u );
                    double _beta = 2. / ( _unorm * _unorm );

                    double _uv_dot = arma::dot( _u, _v );

                    _res = v - _beta * _uv_dot * u;

                    return _res;
                }

                arma::mat mult_HA( const arma::mat &u, const arma::mat &A )
                {
                    assert( u.n_cols == 1 && u.n_rows == A.n_rows );

                    arma::mat _res = arma::zeros<arma::mat>( u.n_rows, A.n_cols );

                    for ( int q = 0; q < A.n_cols; q++ )
                    {
                        _res.col( q ) = mult_Hv( u, A.col( q ) );
                    }

                    return _res;
                }

                arma::mat mult_AH( const arma::mat &u, const arma::mat &A )
                {
                    assert( u.n_cols == 1 && A.n_cols == u.n_rows );

                    arma::mat _res = arma::zeros<arma::mat>( A.n_rows, u.n_rows );

                    for ( int q = 0; q < A.n_rows; q++ )
                    {
                        _res.row( q ) = trans( mult_Hv( u, trans( A.row( q ) ) ) );
                    }

                    return _res;
                }

                ZeroHcalc zero_H( const arma::mat &v )
                {
                    assert( v.n_cols == 1 );

                    ZeroHcalc _res;

                    arma::mat _u;
                    double _sigma = 0.0;

                    arma::vec _v = arma::vectorise( v );
                    double _m = arma::norm( _v, "inf" );

                    _u = v * ( 1. / _m );
                    _sigma = Sign( _u( 0, 0 ) ) * arma::norm( _u, 2 );
                    _u( 0, 0 ) += _sigma;
                    _sigma = - _m * _sigma;

                    _res.u = _u;
                    _res.sigma = _sigma;

                    return _res;
                }

                arma::mat hFromVector( const arma::mat &u )
                {
                	arma::mat _res( u.n_rows, u.n_rows );

                	arma::vec _u = vectorise( u );
                	double _utu = dot( _u, _u );

                	_res = arma::eye<arma::mat>( _res.n_rows, _res.n_cols ) -
                		   2 * u * arma::trans( u ) / _utu;

                	return _res;
                }

                QRcalc factorization( const arma::mat &A )
                {
                    assert( A.n_rows == A.n_cols );

                    int N = A.n_rows;

                    QRcalc _res;

                    arma::mat _Q;
                    arma::mat _R;

                    std::vector<arma::mat> _uu;
                    arma::mat _Ak = A;

                    for ( int q = 0; q < N - 1; q++ )
                    {
                    	// Zero-ize column q of _Ak
                    	ZeroHcalc _res_uk = zero_H( _Ak.submat( q, q, N - 1, q ) );
                    	//cout << "zeroH" << endl;
                    	//cout << _res_uk.u << endl;
                    	_Ak.submat( q, q, N - 1, q ) = arma::zeros<arma::mat>( N - q , 1 );
                    	arma::mat _u = arma::zeros<arma::mat>( N, 1 );
                    	_u.submat( q, 0, N - 1 , 0 ) = _res_uk.u;
                    	_uu.push_back( _u );// Store the Householder related vector
                    	arma::mat _Ak_sub = mult_HA( _res_uk.u, _Ak.submat( q, q, N - 1, N - 1 ) );
                    	_Ak.submat( q, q, N - 1, N - 1 ) = _Ak_sub;
                    	_Ak( q, q ) = _res_uk.sigma;// ?? necessary

                    	//cout << "Ak" << endl;
                    	//cout << _Ak << endl;
                    }

                    _R = _Ak;
                    _Q = arma::eye<arma::mat>( N, N );
                    for ( int q = 0; q < _uu.size(); q++ )
                    {
                    	//cout << "dump vector" << endl;
                    	//cout << _uu[q] << endl;
                    	//cout << "dump matrix" << endl;
                    	//cout << hFromVector( _uu[q] ) << endl;

                    	_Q = _Q * hFromVector( _uu[q] );
                    }

                    _res.Q = _Q;
                    _res.R = _R;

                    return _res;
                }
            }


            namespace gram_schmidt 
            {
            	namespace method
            	{
            		enum _method
            		{
            			CLASSIC,
            			MODIFIED
            		};
            	}

            	vector<arma::vec> gs_orthonormalization( const vector<arma::vec> &vv )
            	{
            		vector<arma::vec> _ee;

            		for ( int q = 0; q < vv.size(); q++ )
            		{
            			arma::vec _u = vv[q];
            			for ( int p = 0; p < q; p++ )
            			{
            				_u -= arma::dot( vv[q], _ee[p] ) * _ee[p];
            			}
            			arma::vec _e = arma::normalise( _u );
            			_ee.push_back( _e );
            		}

            		return _ee;
            	}

            	vector<arma::vec> mgs_orthonormalization( const vector<arma::vec> &vv )
            	{
            		vector<arma::vec> _ee;

            		for ( int q = 0; q < vv.size(); q++ )
            		{
            			arma::vec _u = vv[q];
            			for ( int p = 0; p < q; p++ )
            			{
            				_u -= arma::dot( _u, _ee[p] ) * _ee[p];
            			}
            			arma::vec _e = arma::normalise( _u );
            			_ee.push_back( _e );
            		}

            		return _ee;
            	}            	

            	QRcalc factorization( const arma::mat &A, method::_method gs_method = method::CLASSIC )
            	{
                    assert( A.n_rows == A.n_cols );

                    int N = A.n_rows;

                    QRcalc _res;

                    arma::mat _Q = arma::zeros<arma::mat>( N, N );
                    arma::mat _R = arma::zeros<arma::mat>( N, N );


                    // Generate the orthonormalization of the columns of A
                    vector<arma::vec> _aa;
                    for ( int q = 0; q < N; q++ )
                    {
                    	_aa.push_back( arma::vectorise( A.col( q ) ) );
                    }

                    vector<arma::vec> _ee;
                   	if ( gs_method == method::CLASSIC )
                   	{
                   		_ee = gs_orthonormalization( _aa );
                   	}
                   	else
                   	{
                    	_ee = mgs_orthonormalization( _aa );
                   	}

                    // Gemerate the matrices according to these vects
                   	for ( int col = 0; col < N; col++ )
                   	{
                   		arma::vec _e = _ee[col];
                   		for ( int row = 0; row < N; row++ )
                   		{
                   			_Q( row, col ) = _e( row );
                   			if ( row > col )
                   			{
                   				_R( row, col ) = 0;
                   			}
                   			else
                   			{
	                   			_R( row, col ) = dot( _aa[col], _ee[row] );
                   			}
                   		}
                   	}



                    _res.Q = _Q;
                    _res.R = _R;

                    return _res;
            	}
            }

            namespace givens
            {

            }


        }

    }
}