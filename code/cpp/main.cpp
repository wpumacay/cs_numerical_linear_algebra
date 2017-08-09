
#include <iostream>
#include <armadillo>

#include "HeatPlateSolver.h"

using namespace std;
using namespace arma;


#define LX 1.0
#define LY 1.0

#define DX 0.25
#define DY 0.25

int main()
{

    cout << "initialize system ********" << endl;

    plate2Dsolver::Plate2D plate( LX, LY, DX, DY ); 

    plate.solve();


    return 0;
}
