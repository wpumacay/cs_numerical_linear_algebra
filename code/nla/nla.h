
#pragma once

#include <cassert>
#include <iostream>
#include <armadillo>
#include <vector>
#include <cmath>
#include <algorithm>

#include "nla_lu.h"
#include "nla_qr.h"

using namespace std;
using namespace arma;

namespace NLA
{

	namespace decompositions
	{

		namespace lu
		{
            // Namespace for lu factorizations
		}

		namespace qr
		{
            namespace householder
            {
                // Namespace for the householder functionality ...
                // for the qr factorization
            }


            namespace gram_schmidt
            {
                // Namespace for the gram_schmidt functionality ...
                // for the qr factorization
            }


            namespace givens
            {
                // Namespace for the givens functionality ...
                // for the qr factorization
            }
		}

	}

}


