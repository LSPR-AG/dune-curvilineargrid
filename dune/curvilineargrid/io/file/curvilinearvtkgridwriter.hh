#ifndef DUNE_CURVILINEARVTKGRIDWRITER_HH
#define DUNE_CURVILINEARVTKGRIDWRITER_HH

#include <dune/curvilineargrid/io/file/curvilinearvtkwriter.hh>

/** \brief This is only needed to define the VTK Function at the moment **/
//#include <dune/grid/io/file/vtk/vtkwriter.hh>
//#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

namespace Dune
{



template <class Grid, int mydim>
class VTKScalarFunction
{

public:

	static const int dimension    = Grid::dimension;
	static const int mydimension  = mydim;
	typedef typename Grid::ctype CoordinateType;

	static const bool USE_INTERSECTION = (mydimension == 2);

	// Allow the basis functions to be initialized not only by entities, but also by intersections
	typedef typename Grid::template Codim<dimension - mydimension>::Entity GridEntity;
	typedef typename Grid::LeafGridView::Traits::Intersection  GridIntersection;
	typedef typename std::conditional<USE_INTERSECTION, GridIntersection, GridEntity>::type Entity;

	typedef Dune::FieldVector<CoordinateType, mydimension>  Domain;
	typedef CoordinateType                                  Range;

public:

	VTKScalarFunction() { }

	virtual void init (const Entity & entity) = 0;

	virtual Range evaluate(const Domain & x) const = 0;

	virtual std::string name () const = 0;

};


template <class Grid, int mydim>
class VTKVectorFunction
{
public:

	static const int dimension    = Grid::dimension;
	static const int mydimension  = mydim;
	typedef typename Grid::ctype CoordinateType;

	static const bool USE_INTERSECTION = (mydimension == 2);

	// Allow the basis functions to be initialized not only by entities, but also by intersections
	typedef typename Grid::template Codim<dimension - mydimension>::Entity GridEntity;
	typedef typename Grid::LeafGridView::Traits::Intersection  GridIntersection;
	typedef typename std::conditional<USE_INTERSECTION, GridIntersection, GridEntity>::type Entity;

	typedef Dune::FieldVector<CoordinateType, mydimension>  Domain;
	typedef Dune::FieldVector<CoordinateType, dimension>    Range;

public:

	VTKVectorFunction() { }

	virtual void init (const Entity & entity) = 0;

	virtual Range evaluate(const Domain & x) const = 0;

	virtual std::string name () const = 0;

};







template <class GridType>
class CurvilinearVTKGridWriter
{
    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = 3;
    static const int   EDGE_CODIM     = 2;
    static const int   FACE_CODIM     = 1;
    static const int   ELEMENT_CODIM  = 0;

	typedef Dune::CurvilinearVTKWriter<GridType>  BaseWriter;

public:
	typedef GridType                      Grid;
	static const int dimension = Grid::dimension;

	typedef typename Grid::ctype          CoordinateType;
	typedef typename Grid::LevelGridView  LevelGridView;
	typedef typename Grid::LeafGridView   LeafGridView;

	typedef typename GridType::template Codim<ELEMENT_CODIM>::Entity                         EntityElement;
	typedef typename GridType::template Codim<FACE_CODIM>::Entity                            EntityFace;
	typedef typename GridType::GridBaseType::template Codim<ELEMENT_CODIM>::EntityGeometry   ElementGeometry;
	typedef typename GridType::GridBaseType::template Codim<FACE_CODIM>::EntityGeometry      FaceGeometry;
	typedef typename ElementGeometry::GlobalCoordinate  GlobalCoordinate;

	typedef typename Grid::StructuralType            StructuralType;
	typedef typename Grid::PhysicalTagType           PhysicalTagType;
	typedef typename Grid::InterpolatoryOrderType    InterpolatoryOrderType;

	typedef VTKScalarFunction<Grid, 2>  VTKScalarFunction2D;
	typedef VTKScalarFunction<Grid, 3>  VTKScalarFunction3D;
	typedef VTKVectorFunction<Grid, 2>  VTKVectorFunction2D;
	typedef VTKVectorFunction<Grid, 3>  VTKVectorFunction3D;

	typedef typename std::vector<VTKScalarFunction2D *> VTKScPtrVector2D;
	typedef typename std::vector<VTKScalarFunction3D *> VTKScPtrVector3D;
	typedef typename std::vector<VTKVectorFunction2D *> VTKVecPtrVector2D;
	typedef typename std::vector<VTKVectorFunction3D *> VTKVecPtrVector3D;


	CurvilinearVTKGridWriter(const Grid & grid) :
		vtkScalarFaceFunctionSet_(nullptr),
		vtkScalarElementFunctionSet_(nullptr),
		vtkVectorFaceFunctionSet_(nullptr),
		vtkVectorElementFunctionSet_(nullptr),
		grid_(grid),
		virtualRefinementOrder_(0),
		writePB_(false),
		writeDB_(false),
		writeIB_(false),
		writeGhost_(false),
		writeExplode_(false),
		writePatience_(true)
	{

	}


	// Allow user to fix virtual refiniment level
	void useFixedVirtualRefinement(int newOrder)  { virtualRefinementOrder_ = newOrder; }

	// Allow user to only write interior elements to accelerate writer
	void writeProcessBoundary(bool write)  { writePB_ = write; }
	void writeDomainBoundary(bool write)  { writeDB_ = write; }
	void writeInteriorBoundary(bool write)  { writeIB_ = write; }
	void writeGhost(bool write)						{ writeGhost_ = write; }

	// Allow user to only write interior elements to accelerate writer
	void writePatience(bool patience)  { writePatience_ = patience; }

	// User can choose to use explosion plotting for meshes
	void writeExplode(bool explode)  { writeExplode_ = explode; }



	// Add a scalar or a vector field set to the grid
	// NOTE: IT IS VERY IMPORTANT THAT ALL FIELDS HAVE UNIQUE NAME
	void addFieldSet(VTKScPtrVector2D & field)   { vtkScalarFaceFunctionSet_    = &field; }
	void addFieldSet(VTKScPtrVector3D & field)   { vtkScalarElementFunctionSet_ = &field; }
	void addFieldSet(VTKVecPtrVector2D & field)  { vtkVectorFaceFunctionSet_    = &field; }
	void addFieldSet(VTKVecPtrVector3D & field)  { vtkVectorElementFunctionSet_ = &field; }


	// Write Grid to VTK, including vector field set defined over each element
	void write(std::string path, std::string filenamePrefix)
	{
		// Declare the curvilinear vtk writer
		BaseWriter writer(grid_.comm().rank(), grid_.comm().size());

		// Properties
        bool interpolate = true;
        std::vector<bool>  writeCodim { true, true, false, false };  // Use elements and faces for discretization, and do not use entities of other codimensions

        // Initialize all vector fields, in case there are none on this process
		if (vtkScalarFaceFunctionSet_ != nullptr)  { for (unsigned int i = 0; i < vtkScalarFaceFunctionSet_->size(); i++) { writer.initScalarField( *((*vtkScalarFaceFunctionSet_)[i])); } }
		if (vtkVectorFaceFunctionSet_ != nullptr)  { for (unsigned int i = 0; i < vtkVectorFaceFunctionSet_->size(); i++) { writer.initVectorField( *((*vtkVectorFaceFunctionSet_)[i])); } }
		if (vtkScalarElementFunctionSet_ != nullptr)  { for (unsigned int i = 0; i < vtkScalarElementFunctionSet_->size(); i++) { writer.initScalarField( *((*vtkScalarElementFunctionSet_)[i])); } }
		if (vtkVectorElementFunctionSet_ != nullptr)  { for (unsigned int i = 0; i < vtkVectorElementFunctionSet_->size(); i++) { writer.initVectorField( *((*vtkVectorElementFunctionSet_)[i])); } }

		// Iterate over elements
  		/** \brief Iterate ove all elements of Interior Border partition */
        int elemIterCount = 0;         // Count elements
		LeafGridView leafView = grid_.leafGridView();
		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Writing Elements" );
		for (auto&& elementThis : elements(leafView, Dune::Partitions::all))
		{
			if (writePatience_) {
				LoggingMessage::writePatience("CurvilinearVTKGridWriter: Writing Elements to VTK...", elemIterCount++, grid_.size(ELEMENT_CODIM));
			}

			Dune::GeometryType geomtype   = elementThis.type();
			StructuralType     thisPType  = elementThis.partitionType();

			if (writeGhost_ || (thisPType != Dune::PartitionType::GhostEntity))
			{
				// Constructing a geometry is quite expensive, do it only once
				//EntityGeometry geom = it->geometry();
				ElementGeometry geom                    = grid_.template entityBaseGeometry<ELEMENT_CODIM>(elementThis);
				std::vector<GlobalCoordinate>  nodeSet = geom.vertexSet();

				PhysicalTagType         physicalTag  = grid_.template entityPhysicalTag<ELEMENT_CODIM>(elementThis);
				InterpolatoryOrderType  elementOrder = grid_.template entityInterpolationOrder<ELEMENT_CODIM>(elementThis);
				std::vector<int>        tagSet  { physicalTag, thisPType, grid_.comm().rank() };

				// To accelerate vtk writer, do not overdiscretize low order meshes. Make discretization appropriate for element, unless user specifically asks for fixed order
				int nDiscretizationPoint = (virtualRefinementOrder_ > 0) ? virtualRefinementOrder_ : elementOrder + 4;

				// Write element to VTK
				LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: --Inserting element geometry" );
				writer.template addCurvilinearElement<dimension - ELEMENT_CODIM>(geomtype, nodeSet, tagSet, elementOrder, nDiscretizationPoint, interpolate, writeExplode_, writeCodim);


				if (thisPType != Dune::PartitionType::GhostEntity)
				{
					// It is unexpected that the process would know the field on its own Ghost element. It is most likely that for memory reasons this field
					// will be stored only once, and on a process for which this element is an interior element
					writeScalarField(vtkScalarElementFunctionSet_, elementThis, writer);
					writeVectorField(vtkVectorElementFunctionSet_, elementThis, writer);

					// Write faces - Domain and Process Boundaries
					if (writeDB_ || writePB_ || writeIB_) // Do it only if user wants to see them
					{
						for (auto&& intersection : intersections(leafView, elementThis))
						{
							const EntityFace faceThis = elementThis.template subEntity<FACE_CODIM>(intersection.indexInInside());
							PhysicalTagType faceTag = grid_.template entityPhysicalTag<FACE_CODIM>(faceThis);
							Dune::GeometryType  faceGeomType   = faceThis.type();
							StructuralType thisFacePType  = faceThis.partitionType();

							bool isDB = ((intersection.neighbor() == false) && (intersection.boundary() == true));
							bool isPB = ((intersection.neighbor() == true) && (intersection.outside().partitionType() == Dune::PartitionType::GhostEntity));
							bool isIB = ((!isDB) && (faceTag >= 0));


							if  (( isDB && writeDB_ ) || ( isPB && writePB_ ) || (isIB && writeIB_)) {

											if (isDB)  { thisFacePType = CurvGrid::BOUNDARY_SEGMENT_PARTITION_TYPE; }
								else		if (isIB)   { thisFacePType = CurvGrid::INTERIOR_BOUNDARY_SEGMENT_PARTITION_TYPE; }

								// Constructing a geometry is quite expensive, do it only once
								//EntityGeometry geom = it->geometry();
								FaceGeometry facegeom = grid_.template entityBaseGeometry<FACE_CODIM>(faceThis);
								std::vector<GlobalCoordinate>  faceNodeSet = facegeom.vertexSet();

								std::vector<int>        faceTagSet  { faceTag, thisFacePType, grid_.comm().rank() };

								// Write element to VTK
								LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: --Inserting face geometry" );
								writer.template addCurvilinearElement<dimension - FACE_CODIM>(faceGeomType, faceNodeSet, faceTagSet, elementOrder, nDiscretizationPoint, interpolate, writeExplode_, writeCodim);

								// For now, only write face-based fields for domain boundaries
								if (isDB)  {
									writeScalarField(vtkScalarFaceFunctionSet_, intersection, writer);
									writeVectorField(vtkVectorFaceFunctionSet_, intersection, writer);
								}
							}
						}
					}
				}

			}
		}

		// Write the data to vtk file
		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Writing .vtk file = " + path + filenamePrefix + ".pvtu");

		//writer_.writeVTK(filenamePrefix + ".vtk");
		writer.writeParallelVTU(path, filenamePrefix);

		// With current design it is expected that the VTK Writer deletes all the VTK functions
		deleteField(vtkScalarFaceFunctionSet_);
		deleteField(vtkVectorFaceFunctionSet_);
		deleteField(vtkScalarElementFunctionSet_);
		deleteField(vtkVectorElementFunctionSet_);
	}




private:


	//template <class VTKVecPtrVector, class EntityType>
	template <class VTKFunction, class EntityType>
	void writeScalarField(std::vector<VTKFunction *> * vec, const EntityType & entity, BaseWriter & writer) {
		if (vec != nullptr) {
			for (unsigned int i = 0; i < vec->size(); i++)
			{
				LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Computing element scalar field" );
				(*vec)[i]->init(entity);

				LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Inserting element scalar field" );
				writer.template addScalarField<VTKFunction>(*((*vec)[i]));
			}
		}
	}


	template <class VTKFunction, class EntityType>
	void writeVectorField(std::vector<VTKFunction *> * vec, const EntityType & entity, BaseWriter & writer) {
		if (vec != nullptr) {
			for (unsigned int i = 0; i < vec->size(); i++)
			{
				LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Computing element vector field" );
				(*vec)[i]->init(entity);

				LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Inserting element vector field" );
				writer.template addVectorField<VTKFunction>(*((*vec)[i]));
			}
		}
	}


	template <class VTKVecPtrVector>
	void deleteField(VTKVecPtrVector * vec) const {
		if (vec != nullptr) {
			for (unsigned int i = 0; i < vec->size(); i++)  { delete (*vec)[i]; }
			vec->clear();
		}
	}


	VTKScPtrVector2D  * vtkScalarFaceFunctionSet_;      // Store scalar functions defined over faces
	VTKVecPtrVector2D * vtkVectorFaceFunctionSet_;      // Store vector functions defined over faces
	VTKScPtrVector3D  * vtkScalarElementFunctionSet_;   // Store scalar functions defined over elements
	VTKVecPtrVector3D * vtkVectorElementFunctionSet_;   // Store vector functions defined over elements


	const Grid & grid_;
	int virtualRefinementOrder_;   // User-defined discretization order for element sub-refinement
	bool writePB_;       		// User can choose to write Process Boundary Segments [Default - false]
	bool writeDB_;       		// User can choose to write Domain Boundary Segments [Default - false]
	bool writeIB_;       		// User can choose to write Interior Boundary Segments [Default - false]
	bool writeGhost_;       		// User can choose to write Ghost Elements [Default - false]
	bool writePatience_;           // User can choose to not write patience output for aestetical purposes
	bool writeExplode_;            // User can choose to use explosion plotting for meshes

};




} // Namespace Dune

#endif // DUNE_CURVILINEARVTKGRIDWRITER_HH
