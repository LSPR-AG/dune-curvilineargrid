#ifndef DUNE_GMSH_2_DUNE_MAPPER_HH_
#define DUNE_GMSH_2_DUNE_MAPPER_HH_

namespace Dune
{


// GMSH Triangle Strategy
//   1) Vertices {0,1,2}
//   2) Edges {(01), (12), (20)}
//   3) Recursive internal triangle

// GMSH Tetrahedron Strategy
//   1) Vertices {0,3,1,2}
//   2) Edges {01,12,20,30,32,31}
//   3) Recursive Faces {021, 013, 032, 312}
//   4) Recursive internal tetrahedron

// TODO: Write analytic routine, and use less complicated recursive transformation
// GMSH Tetrahedron Strategy
//   1) Vertices {0,1,2,3}   3->1  1->2  2->3
//   2) Edges {02,23,30,10,13,12}
//   3) Recursive Faces {032, 021, 013, 123}
//   4) Recursive internal tetrahedron

class Gmsh2DuneMapper
{

public:

Gmsh2DuneMapper()
{
    // Initialize triangular point renumberings for GMSH->DUNE convention
    triangularInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 1, 2} );
    triangularInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 3, 1, 5, 4, 2} );
    triangularInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 3, 4, 1, 8, 9, 5, 7, 6, 2} );
    triangularInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 3, 4, 5, 1, 11, 12, 13, 6, 10, 14, 7, 9, 8, 2} );
    triangularInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 3, 4, 5, 6, 1, 14, 15, 18, 16, 7, 13, 20, 19, 8, 12, 17, 9, 11, 10, 2} );

    // Initialize tetrahedral point renumberings for GMSH->DUNE convention
    tetrahedralInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 3, 1, 2} );
    tetrahedralInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 7, 3, 4, 9, 1, 6, 8, 5, 2} );
    tetrahedralInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 11, 10, 3, 4, 17, 14, 5, 15, 1, 9, 18, 12, 16, 19, 6, 8, 13, 7, 2} );
    tetrahedralInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 15, 14, 13, 3, 4, 25, 27, 19, 5, 26, 20, 6, 21, 1, 12, 28, 29, 16, 22, 34, 31, 24, 32, 7, 11, 30, 17, 23, 33, 8, 10, 18, 9, 2} );
    tetrahedralInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 19, 18, 17, 16, 3, 4, 34, 39, 36, 24, 5, 37, 38, 25, 6, 35, 26, 7, 27, 1, 15, 40, 43, 41, 20, 28, 52, 55, 46, 33, 53, 49, 30, 47, 8, 14, 45, 44, 21, 31, 54, 51, 32, 50, 9, 13, 42, 22, 29, 48, 10, 12, 23, 11, 2} );
}


// ***********************************************************************
// GMSH Convention Methods
// ***********************************************************************

// Constructs a DUNE geometry type based on GMSH element index
GeometryType geometryType(int gmshIndex)
{
    GeometryType rez;
    int gi = gmshIndex;

         if (gi == 15)                                                                                                     { rez.makeVertex(); }
    else if ((gi == 1) || (gi == 8) || (gi == 26) || (gi == 27) || (gi == 28))                                             { rez.makeLine(); }
    else if ((gi == 2) || (gi == 9) || (gi == 20) || (gi == 21) || (gi == 22) || (gi == 23) || (gi == 24) || (gi == 25))   { rez.makeTriangle(); }
    else if ((gi == 3) || (gi == 10) || (gi == 16))                                                                        { rez.makeQuadrilateral(); }
    else if ((gi == 4) || (gi == 11) || (gi == 29) || (gi == 30) || (gi == 31))                                            { rez.makeTetrahedron(); }
    else if ((gi == 5) || (gi == 12) || (gi == 17) || (gi == 92) || (gi == 93))                                            { rez.makeHexahedron(); }
    else if ((gi == 6) || (gi == 13) || (gi == 18))                                                                        { rez.makePrism(); }
    else if ((gi == 7) || (gi == 14) || (gi == 19))                                                                        { rez.makePyramid(); }
    else  { DUNE_THROW(Dune::IOError, "Unexpected geometry type");  }

/*             // Note that hexahedron GeometryType is missing
         switch (thisElmName)
         {
         case GMSH_EDGE           : factory.insertElement(GeometryType(GeometryType::simplex,dim_), corners);   break;
         case GMSH_TRIANGLE       : factory.insertElement(GeometryType(GeometryType::simplex,dim_), corners);   break;
         case GMSH_QUADRANGLE     : factory.insertElement(GeometryType(GeometryType::cube,dim_),    corners);   break;
         case GMSH_TETRAHEDRON    : factory.insertElement(GeometryType(GeometryType::simplex,dim_), corners);   break;
         case GMSH_HEXAHEDRON     : factory.insertElement(GeometryType(GeometryType::simplex,dim_), corners);   break;
         case GMSH_PRISM          : factory.insertElement(GeometryType(GeometryType::prism,dim_),   corners);   break;
         case GMSH_PYRAMID        : factory.insertElement(GeometryType(GeometryType::pyramid,dim_), corners);   break;
         }*/

    return rez;
}


// Returns the type name of the element given its GMSH_index
int elementOrder(int gmshIndex)
{
    // Hexahedra of high dimension have funny index
    if (gmshIndex == 92) { return 3; }
    if (gmshIndex == 93) { return 4; }

    // Array copy-pasted from GMSH Brute-Force because it does not seem to have any pattern :)
    const int elemOrder[32]          = {1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 3, 4, 5, 3, 4, 5};

    return elemOrder[gmsh2DuneIndex(gmshIndex)];
}


// Tells which GMSH indices point to entities of incomplete polynomial order
bool hasIncompleteOrder(int gmshIndex)
{
    // Both high-order hexahedrons are complete
    if ((gmshIndex == 92) || (gmshIndex == 93)) { return 0; }

    // Array copy-pasted from GMSH Brute-Force because it does not seem to have any pattern :)
    const bool elemIncomplete[32]          = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0};

    return elemIncomplete[gmsh2DuneIndex(gmshIndex)];
}


// Returns the number of degrees of freedom of the element given its GMSH_index
// note: This info can not simply be obtained from referenceElement, because some of the elements in GMSH have incomplete order, so less DoF than expected
int dofNumber(int gmshIndex)
{
    // Array copy-pasted from GMSH Brute-Force because it does not seem to have any pattern :)
    const int nDofs[32]              = {2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, 15, 15, 21, 4, 5, 6, 20, 35, 56};

    return nDofs[gmsh2DuneIndex(gmshIndex)];
}


// Returns the total number of DoF associated with all subentities of a given dimension for this element, subtracting the ones that come from the corners
int subentityExtraDofNumber(int gmshIndex, int dim)
{
    const int nDofsExtraEdge[32] = {0, 0, 0, 0, 0, 0, 0, 1, 3, 4, 6, 12, 9, 8, 0, 4, 12, 9, 8, 6, 6, 9, 9, 12, 12, 2, 3, 4, 12, 18, 24};
    const int nDofsExtraFace[32] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 6, 3, 1, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 6, 0, 0, 0, 4, 12, 24};
    const int nDofsExtraElem[32] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4};

    switch(dim)
    {
    case 1: return nDofsExtraEdge[gmsh2DuneIndex(gmshIndex)];  break;
    case 2: return nDofsExtraFace[gmsh2DuneIndex(gmshIndex)];  break;
    case 3: return nDofsExtraElem[gmsh2DuneIndex(gmshIndex)];  break;
    }

    return -1;
}


// correct differences between gmsh and Dune in the local vertex numbering
// [FIXME] THIS METHOD DOES NOT WORK WITH INCOMPLETE ORDER GMSH ELEMENTS AT THE MOMENT
void gmsh2DuneElementDofNumbering(GeometryType gt, int thisElmOrder, std::vector<int> &elementDofSet) {
    int thisElmDofNo = elementDofSet.size();
    std::vector<int> tmp;

    if (!gt.isSimplex() || gt.dim() < 2 || gt.dim() > 3)  {
    	DUNE_THROW(Dune::IOError, "CURVILINEAR_GMSH_READER: gmsh2duneRenumbering() only implemented for Simplex 2D and 3D geometries at the moment");
    }

    if (gt.isTriangle())
    {
        for (int i = 0; i < thisElmDofNo; i++) { tmp.push_back(elementDofSet[triangularInterpolatoryVertexGmsh2DuneMap[thisElmOrder - 1][i]]); }
        for (int i = 0; i < thisElmDofNo; i++) { elementDofSet[i] = tmp[i]; }
    }
    else if (gt.isTetrahedron())
    {
        for (int i = 0; i < thisElmDofNo; i++) { tmp.push_back(elementDofSet[tetrahedralInterpolatoryVertexGmsh2DuneMap[thisElmOrder - 1][i]]); }
        for (int i = 0; i < thisElmDofNo; i++) { elementDofSet[i] = tmp[i]; }
    }
}


// In GMSH the global vertex index starts at 1, in Dune it starts at 0, therefore correction
int  gmsh2DuneIndex (int gmshIndex) { return gmshIndex - 1; }








private:
    // A map from GMSH -> Dune for indexing interpolatory points
    std::vector< std::vector< int > > triangularInterpolatoryVertexGmsh2DuneMap;
    std::vector< std::vector< int > > tetrahedralInterpolatoryVertexGmsh2DuneMap;
};


} // namespace Dune

#endif // DUNE_GMSH_2_DUNE_MAPPER_HH_
