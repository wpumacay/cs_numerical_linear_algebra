
#include <iostream>
#include "nla/nla.h"


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
    cout << "C: " << endl;
    cout << C << endl;
    cout << "H: " << endl;
    cout << _h << endl;
    cout << "H-arma: " << endl;
    cout << arma::chol( C, "lower" ) << endl;


    return 0;
}