// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>

#include <dune/curvilineargrid/utility/paralleldatawriter.hh>





/** \brief
 *
 * The following tutorial demonstrates the use of intrinsic global index provided by curvilinear grid.
 * In particular, a special parameter of the grid and GMSH reader allows to construct the global index for elements (codim 0) by reusing the GMSH intrinsic index.
 * The advantage of such global index is that same element has the same global index, regardless of the number of cores the mesh is constructed on, which is very useful for debugging parallel FEM codes.
 *
 * The below example uses a small utility called ParallelDataWriter to write the volumes of all elements of the grid in a single file, sorted by their global index.
 *
 * The user could verify by himself (e.g. using Matlab) that this code produces exactly the same output file regardless of the number of processes it is run on.
 */


int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;

	// Create Grid
	typedef Dune::CurvilinearGrid<ctype, dim, isCached, Dune::LoggingMessage> GridType;
	GridType * grid = createGrid<GridType>(mpihelper);

	const bool isCached = true;
	const int ELEMENT_CODIM = 0;  // Codimension of element in 3D

	// Reporting vector
	typedef int     IndexType;
	typedef double  DataType;
	typedef Dune::ParallelDataWriter<GridType, IndexType, DataType>  PDWVector;

	std::vector<IndexType>  indexVec;
	std::vector<int>        sizeVec;
	std::vector<DataType>   dataVec;

	int elemInd = 0;
	int nElem = grid->numInternal(ELEMENT_CODIM);

	typename GridType::LeafGridView leafView = grid->leafGridView();
	for (auto&& elementThis : elements(leafView, Dune::Partitions::interiorBorder))
	{
		Dune::LoggingMessage::writePatience("iterating over elements ", elemInd++, nElem);

		int globalIndex = grid->template entityGlobalIndex<ELEMENT_CODIM>(elementThis);
		indexVec.push_back(globalIndex);
		sizeVec.push_back(1);
		dataVec.push_back(elementThis.geometry().volume());
	}

	PDWVector::writeParallelData2File("volumevector_" + std::to_string(mpihelper.size()) + ".txt", indexVec, sizeVec, dataVec, *grid);

	typedef Dune::LoggingTimer<Dune::LoggingMessage>                 LoggingTimerDev;
	LoggingTimerDev::getInstance().reportParallel(mpihelper);

    // Delete the grid
    delete grid;


    return 0;
}
