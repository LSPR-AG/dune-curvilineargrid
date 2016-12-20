// Autotool generated header.
#include <config.h>

// Stl headers.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>



using namespace Dune;

using namespace Dune::CurvGrid;


int main(int argc, char** argv)
{

	std::vector<int> A {6,5,3,6,3,6,1,34,55,6,45,2,1,4};
	std::vector<int> Ac = A;
	VectorHelper::compactify(Ac);
	std::vector<int> B {3,4,5,6,7,8};

	std::vector<int> AcplusB = VectorHelper::sortedSetUnion(Ac, B);
	std::vector<int> AcminusB = VectorHelper::sortedSetComplement(Ac, B);
	std::vector<int> AcandB = VectorHelper::sortedSetIntersection(Ac, B);

	std::cout << " start with vector A=" + VectorHelper::vector2string(A) << std::endl;
	std::cout << " its compactification is Ac =" << VectorHelper::vector2string(Ac) << std::endl;
	std::cout << " also consider vector B =" << VectorHelper::vector2string(B) << std::endl;
	std::cout << " the union      (Ac + B) = " << VectorHelper::vector2string(AcplusB) << std::endl;
	std::cout << " the complement (Ac - B) = " << VectorHelper::vector2string(AcminusB) << std::endl;
	std::cout << " the intersect  (Ac & B) = " << VectorHelper::vector2string(AcandB) << std::endl;

    /** \brief leave program peacefully */
    return(0);
}
