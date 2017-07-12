
#include <iostream>
#include "nla/nla.h"

using namespace std;
using namespace NLA;

int main()
{
    cout << "testing householder functionality ***********" << endl;

    cout << "zero_H *******************" << endl;


    arma::mat x;
    x << -0.2071 << arma::endr
      << -1.2071 << arma::endr;

    decompositions::qr::ZeroHcalc _res_1 = decompositions::qr::householder::zero_H( x );

    cout << "x" << endl;
    cout << x << endl;

    cout << "u" << endl;
    cout << _res_1.u << endl;
    cout << "sigma" << endl;
    cout << _res_1.sigma << endl;

    cout << "**************************" << endl;

    cout << "qr-householder ***********" << endl;

    arma::mat A = { { 0,  1,  1 },
                    { 1,  2,  3 },
                    { 1,  1,  1 } };

    decompositions::qr::QRcalc _res_2 = decompositions::qr::householder::factorization( A );

    cout << "A" << endl;
    cout << A << endl;

    cout << "Q" << endl;
    cout << _res_2.Q << endl;
    cout << "R" << endl;
    cout << _res_2.R << endl;
    cout << "QR" << endl;
    cout << _res_2.Q * _res_2.R << endl;
    cout << "QQ^{T}" << endl;
    cout << _res_2.Q * trans( _res_2.Q ) << endl;

    cout << "|| I - QQ^{T} ||" << endl << endl;
    cout << arma::norm( arma::eye<arma::mat>( 3, 3 ) - _res_2.Q * trans( _res_2.Q ) ) << endl;

    cout << "**************************" << endl;

    cout << "qr-gram_schmidt classic **" << endl;

    decompositions::qr::QRcalc _res_3 = decompositions::qr::gram_schmidt::factorization( A );

    cout << "A" << endl;
    cout << A << endl;

    cout << "Q" << endl;
    cout << _res_3.Q << endl;
    cout << "R" << endl;
    cout << _res_3.R << endl;
    cout << "QR" << endl;
    cout << _res_3.Q * _res_3.R << endl;
    cout << "QQ^{T}" << endl;
    cout << _res_3.Q * trans( _res_3.Q ) << endl;
    cout << "|| I - QQ^{T} ||" << endl << endl;
    cout << arma::norm( arma::eye<arma::mat>( 3, 3 ) - _res_3.Q * trans( _res_3.Q ) ) << endl;

    cout << "**************************" << endl;

    cout << "qr-gram_schmidt modifed **" << endl;

    decompositions::qr::QRcalc _res_4 = 
            decompositions::qr::gram_schmidt::factorization( A, 
                                                             decompositions::qr::gram_schmidt::method::MODIFIED );

    cout << "A" << endl;
    cout << A << endl;

    cout << "Q" << endl;
    cout << _res_4.Q << endl;
    cout << "R" << endl;
    cout << _res_4.R << endl;
    cout << "QR" << endl;
    cout << _res_4.Q * _res_4.R << endl;
    cout << "QQ^{T}" << endl;
    cout << _res_4.Q * trans( _res_4.Q ) << endl;
    cout << "|| I - QQ^{T} ||" << endl << endl;
    cout << arma::norm( arma::eye<arma::mat>( 3, 3 ) - _res_4.Q * trans( _res_4.Q ) ) << endl;

    cout << "**************************" << endl;

    cout << "qr hessenberg test mat ***" << endl;

    int N = 10;
    arma::mat H = arma::zeros<arma::mat>( N, N );
    for ( int col = 0; col < N; col++ )
    {
        for ( int row = 0; row < N; row++ )
        {
            H( row, col ) = 1. / ( row + col + 1 );
        }        
    }

    decompositions::qr::QRcalc _res_5 = decompositions::qr::gram_schmidt::factorization( H );
    //cout << "QQ^{T}" << endl;
    //cout << _res_5.Q * trans( _res_5.Q ) << endl;
    cout << "|| I - QQ^{T} ||" << endl << endl;
    cout << arma::norm( arma::eye<arma::mat>( N, N ) - _res_5.Q * trans( _res_5.Q ) ) << endl;

    decompositions::qr::QRcalc _res_6 = 
            decompositions::qr::gram_schmidt::factorization( H,
                                                             decompositions::qr::gram_schmidt::method::MODIFIED );
    //cout << "QQ^{T}" << endl;
    //cout << _res_6.Q * trans( _res_6.Q ) << endl;
    cout << "|| I - QQ^{T} ||" << endl << endl;
    cout << arma::norm( arma::eye<arma::mat>( N, N ) - _res_6.Q * trans( _res_6.Q ) ) << endl;

    decompositions::qr::QRcalc _res_7 = decompositions::qr::householder::factorization( H );
    //cout << "QQ^{T}" << endl;
    //cout << _res_6.Q * trans( _res_6.Q ) << endl;
    cout << "|| I - QQ^{T} ||" << endl << endl;
    cout << arma::norm( arma::eye<arma::mat>( N, N ) - _res_7.Q * trans( _res_7.Q ) ) << endl;

    cout << "**************************" << endl;

    cout << "*********************************************" << endl;

    return 0;
}