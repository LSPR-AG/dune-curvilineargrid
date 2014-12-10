/***************************************************************************
                          tetmesh.h  -  description
                             -------------------
    begin                : Mon Dec 15 2003
    copyright            : (C) 2003 by Roman Geus
    email                : roman.geus@psi.ch
    edited by            : Hua Guo, Sep 2010
***************************************************************************/

/***************************************************************************
                          curvilineargridbase.hh
                             -------------------
    begin                : Tue Nov 25 2014
    copyright            : (C) 2014 by Aleksejs Fomins, LSPR AG
    description          : Upgraded the mesh to curvilinear grid
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DUNE_CURVILINEARGRIDDIAGNOSTIC_HH
#define DUNE_CURVILINEARGRIDDIAGNOSTIC_HH

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/curvilinearelementinterpolator.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>



namespace Dune {

template <class ct>
class CurvilinearGridDiagnostic
{
private:
	Dune::CurvilinearGridBase<ct> & gridbase_;

public:
	CurvilinearGridDiagnostic(Dune::CurvilinearGridBase<ct> & gridbase) : gridbase_(gridbase)
	{

	}





	// Tests
	// *************************************************************8


	// Writes the mesh to VTK, including the additional constructions made by the mesh generator
	void vtkWriteMesh (
		bool withDomainBoundaries,
		bool withProcessBoundaries,
		bool withGhostElements
	) {}


	// Computes the ratio of shortest vs longest edge for each element
	std::vector<double> linearElementQuality() {}

	std::vector<double> linearElementVolume() {}

	std::vector<double> curvilinearElementVolume() {}

	// Computes the ratio of curvilinear vs linear volume for each element
	std::vector<double> curvilinearElementVolumeRatio() {}

	double processBoundarySurfaceArea() {}

	double domainBoundarySurfaceArea() {}


};


} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDDIAGNOSTIC_HH
