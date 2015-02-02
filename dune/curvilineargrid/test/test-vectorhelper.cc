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


int main(int argc, char** argv)
{

	std::vector<int> A {6,5,3,6,3,6,1,34,55,6,45,2,1,4};
	std::vector<int> Ac = A;
	Dune::VectorHelper::compactify(Ac);
	std::vector<int> B {3,4,5,6,7,8};

	std::vector<int> AcplusB = Dune::VectorHelper::sortedSetUnion(Ac, B);
	std::vector<int> AcminusB = Dune::VectorHelper::sortedSetComplement(Ac, B);
	std::vector<int> AcandB = Dune::VectorHelper::sortedSetIntersection(Ac, B);

	std::cout << " start with vector A=" + Dune::VectorHelper::vector2string(A) << std::endl;
	std::cout << " its compactification is Ac =" << Dune::VectorHelper::vector2string(Ac) << std::endl;
	std::cout << " also consider vector B =" << Dune::VectorHelper::vector2string(B) << std::endl;
	std::cout << " the union      (Ac + B) = " << Dune::VectorHelper::vector2string(AcplusB) << std::endl;
	std::cout << " the complement (Ac - B) = " << Dune::VectorHelper::vector2string(AcminusB) << std::endl;
	std::cout << " the intersect  (Ac & B) = " << Dune::VectorHelper::vector2string(AcandB) << std::endl;

    /** \brief leave program peacefully */
    return(0);
}
