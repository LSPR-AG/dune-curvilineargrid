#ifndef DUNE_CURVILINEARVTKGRIDWRITER_HH
#define DUNE_CURVILINEARVTKGRIDWRITER_HH

#include <dune/curvilineargrid/io/file/curvilinearvtkwriter.hh>

/** \brief This is only needed to define the VTK Function at the moment **/
//#include <dune/grid/io/file/vtk/vtkwriter.hh>
//#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

namespace Dune
{

namespace CurvGrid {


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
    // Dimensions of entity types for better code readability
    static const int   DIM0D   = 0;
    static const int   DIM1D   = 1;
    static const int   DIM2D   = 2;
    static const int   DIM3D   = 3;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = 3;
    static const int   EDGE_CODIM     = 2;
    static const int   FACE_CODIM     = 1;
    static const int   ELEMENT_CODIM  = 0;

	typedef CurvilinearVTKWriter<GridType>  BaseWriter;

public:
	typedef GridType                      Grid;
	static const int dimension = Grid::dimension;

	typedef typename Grid::ctype          CoordinateType;
	typedef typename Grid::LevelGridView  LevelGridView;
	typedef typename Grid::LeafGridView   LeafGridView;

	typedef typename GridType::GlobalIndexType  GlobalIndexType;
	typedef typename GridType::template Codim<ELEMENT_CODIM>::Entity                         EntityElement;
	typedef typename GridType::template Codim<FACE_CODIM>::Entity                            EntityFace;
	typedef typename GridType::LeafGridView::Traits::Intersection  Intersection;
	typedef typename GridType::GridBaseType::GridEntity::template Codim<ELEMENT_CODIM>::EntityGeometry   ElementGeometry;
	typedef typename GridType::GridBaseType::GridEntity::template Codim<FACE_CODIM>::EntityGeometry      FaceGeometry;
	typedef typename ElementGeometry::GlobalCoordinate  GlobalCoordinate;
	typedef std::vector<GlobalCoordinate> GlobalVector;


	typedef Dune::ReferenceElements<CoordinateType, dimension-1> ReferenceElements2d;

	typedef typename Grid::StructuralType            StructuralType;
	typedef typename Grid::PhysicalTagType           PhysicalTagType;
	typedef typename Grid::InterpolatoryOrderType    InterpolatoryOrderType;

	typedef VTKScalarFunction<Grid, DIM2D>  VTKScalarFunction2D;
	typedef VTKScalarFunction<Grid, DIM3D>  VTKScalarFunction3D;
	typedef VTKVectorFunction<Grid, DIM2D>  VTKVectorFunction2D;
	typedef VTKVectorFunction<Grid, DIM3D>  VTKVectorFunction3D;

	typedef typename std::vector<VTKScalarFunction2D *> VTKScPtrVector2D;
	typedef typename std::vector<VTKScalarFunction3D *> VTKScPtrVector3D;
	typedef typename std::vector<VTKVectorFunction2D *> VTKVecPtrVector2D;
	typedef typename std::vector<VTKVectorFunction3D *> VTKVecPtrVector3D;


	CurvilinearVTKGridWriter(const Grid & grid) :
		grid_(grid),
		vtkScalarFaceFunctionSet_(nullptr),
		vtkScalarElementFunctionSet_(nullptr),
		vtkVectorFaceFunctionSet_(nullptr),
		vtkVectorElementFunctionSet_(nullptr),
		virtualRefinementOrder_(0),
		writePB_(false),
		writeDB_(false),
		writeIB_(false),
		writePeriodic_(false),
		writePeriodicBind_(false),
		writeGhost_(false),
		writeInterpolate_(true),
		writeExplode_(false),
		writePatience_(true),
		vtkDataFormat_(VTU_DATA_FORMAT_ASCII)
	{

	}


	// Allow user to fix virtual refiniment level
	void useFixedVirtualRefinement(int newOrder)  { virtualRefinementOrder_ = newOrder; }

	// Allow user to only write interior elements to accelerate writer
	void writeProcessBoundary(bool write)  { writePB_ = write; }
	void writeDomainBoundary(bool write)  { writeDB_ = write; }
	void writeInteriorBoundary(bool write)  { writeIB_ = write; }
	void writePeriodicBoundary(bool write)  { writePeriodic_ = write; }
	void writeGhost(bool write)						{ writeGhost_ = write; }
	void writePeriodicBind(bool write)  { writePeriodicBind_ = write; }

	// Allow user to only write interior elements to accelerate writer
	void writePatience(bool patience)  { writePatience_ = patience; }

	// User can choose to use explosion plotting for meshes
	void writeInterpolate(bool interpolate)  { writeInterpolate_ = interpolate; }
	void writeExplode(bool explode)  { writeExplode_ = explode; }

	// User set VTK output format
	void setVTKDataFormat(std::string format) { vtkDataFormat_ = format; }

	// Add a scalar or a vector field set to the grid
	// NOTE: IT IS VERY IMPORTANT THAT ALL FIELDS HAVE A UNIQUE NAME
	void addFieldSet(VTKScPtrVector2D & field)   { vtkScalarFaceFunctionSet_    = &field; }
	void addFieldSet(VTKScPtrVector3D & field)   { vtkScalarElementFunctionSet_ = &field; }
	void addFieldSet(VTKVecPtrVector2D & field)  { vtkVectorFaceFunctionSet_    = &field; }
	void addFieldSet(VTKVecPtrVector3D & field)  { vtkVectorElementFunctionSet_ = &field; }


	// Write Grid to VTK, including vector field set defined over each element
	void write(std::string path, std::string filenamePrefix)
	{
		// Construct base writer
		// Important note: Basis writer stores a lot of data. Current paradigm is to construct a new writer for each file, and destroy the writer after writing is done.
		BaseWriter baseWriter(grid_.comm().rank(), grid_.comm().size());


		// Determine Grid Capabilities
		bool isGridParallel = grid_.comm().size() > 1;
		bool hasGridGhost = grid_.ghostSize(ELEMENT_CODIM) > 0;

        // Determine if any of the faces will be written
        bool writeAnyFace = writeDB_ || writePB_ || writeIB_ || writePeriodic_;

        // Initialize all vector fields, in case there are none on this process
		if (vtkScalarFaceFunctionSet_ != nullptr)  { for ( const auto & scalarFaceFunctor : *vtkScalarFaceFunctionSet_) { baseWriter.initScalarField(*scalarFaceFunctor); } }
		if (vtkVectorFaceFunctionSet_ != nullptr)  { for ( const auto & vectorFaceFunctor : *vtkVectorFaceFunctionSet_) { baseWriter.initVectorField(*vectorFaceFunctor); } }
		if (vtkScalarElementFunctionSet_ != nullptr)  { for ( const auto & scalarElementFunctor : *vtkScalarElementFunctionSet_) { baseWriter.initScalarField(*scalarElementFunctor); } }
		if (vtkVectorElementFunctionSet_ != nullptr)  { for ( const auto & vectorElementFunctor : *vtkVectorElementFunctionSet_) { baseWriter.initVectorField(*vectorElementFunctor); } }

		// Iterate over elements
  		/** \brief Iterate over all elements of Interior Border partition */
        int elemIterCount = 0;         // Count elements
		LeafGridView leafView = grid_.leafGridView();
		LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Writing Elements" );
		for (auto&& interiorElem : elements(leafView, Dune::Partitions::interiorBorder))
		{
			if (writePatience_) {
				LoggingMessage::writePatience("CurvilinearVTKGridWriter: Writing Elements to VTK...", elemIterCount++, grid_.size(ELEMENT_CODIM));
			}

			// Write this interior element
			Dune::PartitionType interiorElemType = interiorElem.partitionType();
			GlobalIndexType interiorElemGInd = grid_.template entityGlobalIndex<ELEMENT_CODIM>(interiorElem);
			bool interiorElemWriteField = true;

			writeEntity(baseWriter, interiorElem, interiorElem, interiorElemType, interiorElemWriteField);

			// Write faces and ghost elements
			if (writeAnyFace || writeGhost_)
			{
				for (auto&& intersection : intersections(leafView, interiorElem))
				{
					// Get geometric parameters
					const EntityFace intersectionFace = interiorElem.template subEntity<FACE_CODIM>(intersection.indexInInside());
					PhysicalTagType intersectionTag = grid_.template entityPhysicalTag<FACE_CODIM>(intersectionFace);
					GlobalCoordinate intersectionOuterNormal = intersection.unitOuterNormal( ReferenceElements2d::simplex().position( 0, 0 ) );

					// NOTE: CurvGridVTKWriter introduces extended partition type format to better visualise special boundary segments
					int intersectionEffectivePartitionType  = intersectionFace.partitionType();

					// Determine if the face will be written, based on its type
					// Note that periodic boundaries are also domain boundaries
					bool writeFace = false;
					bool isPeriodic = ((intersection.neighbor() == true) && (intersection.boundary() == true));
					bool isDB = ((intersection.neighbor() == false) && (intersection.boundary() == true)) || isPeriodic;
					bool isPB = isGridParallel
							? ((intersection.neighbor() == true) && (intersection.boundary() == false) && (intersection.outside().partitionType() == Dune::PartitionType::GhostEntity))
							: false; // There are no PB in serial case
					bool isIB = ((!isDB) && (intersectionTag >= 0));

//					if (isDB == isPB) {
//						std::stringstream logstr;
//						logstr << "CurvVTKGridWriter: Unexpected intersection setup"
//								<< " hasNB=" << intersection.neighbor()
//								<< " isB=" << intersection.boundary()
//								<< " OutPT=" << intersection.outside().partitionType();
//						DUNE_THROW(Dune::IOError, logstr.str());
//					}

					assert(!(isDB && isPB));				// Can't be DB and PB at the same time

					if				( isPeriodic && writePeriodic_)	{ writeFace = true;  intersectionEffectivePartitionType = PERIODIC_BOUNDARY_PARTITION_TYPE; }
					else if		( isDB && writeDB_)					{ writeFace = true;  intersectionEffectivePartitionType = BOUNDARY_SEGMENT_PARTITION_TYPE; }
					else if		( isIB && writeIB_)						{ writeFace = true;  intersectionEffectivePartitionType = INTERIOR_BOUNDARY_SEGMENT_PARTITION_TYPE; }
					else if		( isPB && writePB_)						{
						// Do not write PB twice, it is wasteful. Determine which process writes it using unit normal at center
						// [TODO] Normal orientation is ugly. Can determine owner natively, if can extract easily neighbor rank
						writeFace = (intersectionOuterNormal[0] > 0) ||
								((intersectionOuterNormal[0] == 0) && (intersectionOuterNormal[1] > 0)) ||
								((intersectionOuterNormal[0] == 0) && (intersectionOuterNormal[1] == 0) && (intersectionOuterNormal[2] > 0));
					}

					// Write face
					// Note: For now, only write face-based fields for domain boundaries
					if  (writeFace) { writeEntity(baseWriter, intersectionFace, intersection, intersectionEffectivePartitionType, isDB); }

					// Write ghost element
					if (writeGhost_ && (isPB || isPeriodic)) {
						EntityElement ghostElement = intersection.outside();
						Dune::PartitionType ghostType = ghostElement.partitionType();
						GlobalIndexType ghostGInd = grid_.template entityGlobalIndex<ELEMENT_CODIM>(ghostElement);

						// Note: It is unexpected that the process would know the field on its own Ghost element.
						// It is most likely that for memory reasons that this field will be stored only once, and on a process for which this element is an interior element
						bool ghostElementWriteField = false;


						if (isPB) {
							assert(ghostType == Dune::PartitionType::GhostEntity);  // In periodic case, periodic ghost can be local interior element. In this case, it is not marked as a ghost
							writeEntity(baseWriter, ghostElement, ghostElement, ghostType, ghostElementWriteField);
						} else {
							writeEntity(baseWriter, ghostElement, ghostElement, PERIODIC_GHOST_PARTITION_TYPE, ghostElementWriteField);

							////////////////////////////////////////////////////////////////////////////////////
							// Test:
							// Write edge connecting CoM of the periodic neighbor faces, as subentities of neighbor elements
							////////////////////////////////////////////////////////////////////////////////////
							if (writePeriodicBind_) {
								GlobalCoordinate interiorCoM = intersection.geometry().center();
								GlobalCoordinate ghostCoM = ghostElement.template subEntity<FACE_CODIM>(intersection.indexInOutside()).geometry().center();
								GlobalVector periodicBindEdge {interiorCoM, ghostCoM};

								Dune::GeometryType edgeType = Dune::GeometryTypes::line;

								// Color the connecting edge by the direction
								int periodicPhysicalTag = -1;
								if (fabs(fabs(intersectionOuterNormal[0]) - 1.0) < NUMERICS_RELATIVE_TOLERANCE) { periodicPhysicalTag = 0;}
								if (fabs(fabs(intersectionOuterNormal[1]) - 1.0) < NUMERICS_RELATIVE_TOLERANCE) { periodicPhysicalTag = 1;}
								if (fabs(fabs(intersectionOuterNormal[2]) - 1.0) < NUMERICS_RELATIVE_TOLERANCE) { periodicPhysicalTag = 2;}

								assert(periodicPhysicalTag != -1);


								// Note: No physical tag for fake edges
								std::vector<int> tagSet {periodicPhysicalTag, PERIODIC_GHOST_BIND_EDGE_TYPE, grid_.comm().rank() };
								int nDiscretizationPoint = 2; // Do not discretize
								InterpolatoryOrderType curvOrder = 1; // Linear

								// Write element to VTK
								std::stringstream logstr;
								logstr << "CurvilinearVTKGridWriter: --Inserting periodic binding edge ";
								LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, logstr.str());

								baseWriter.template addCurvilinearElement<DIM1D, DIM1D>(edgeType, periodicBindEdge, tagSet, curvOrder, nDiscretizationPoint, writeInterpolate_, writeExplode_);
							}
						}
					}
				}
			}
		}

		// Write the data to vtk file
		LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Writing .vtk file = " + path + filenamePrefix + ".pvtu");

		//writer_.writeVTK(filenamePrefix + ".vtk");
		baseWriter.writeParallelVTU(path, filenamePrefix, vtkDataFormat_);

		// With current design it is expected that the VTK Writer deletes all the VTK functions
		deleteField(vtkScalarFaceFunctionSet_);
		deleteField(vtkVectorFaceFunctionSet_);
		deleteField(vtkScalarElementFunctionSet_);
		deleteField(vtkVectorElementFunctionSet_);
	}




protected:

	template <class GridEntityType, class FieldEntityType>
	void writeEntity(BaseWriter & baseWriter, const GridEntityType & gridEntity, const FieldEntityType & fieldEntity, int ptEff, bool withField) {
		const int codim = GridEntityType::codimension;
		const int mydim = dimension - codim;
		Dune::GeometryType gt   = gridEntity.type();

		// Constructing a geometry is quite expensive, do it only once
		typedef typename GridType::GridBaseType::GridEntity::template Codim<codim>::EntityGeometry   ThisGeometry;
		ThisGeometry geom = grid_.template entityBaseGeometry<codim>(gridEntity);
		GlobalVector nodeSet = geom.vertexSet();

		PhysicalTagType physicalTag  = grid_.template entityPhysicalTag<codim>(gridEntity);
		InterpolatoryOrderType curvOrder = grid_.template entityInterpolationOrder<codim>(gridEntity);
		std::vector<int> tagSet { physicalTag, ptEff, grid_.comm().rank() };

		// To accelerate vtk writer, do not overdiscretize low order meshes. Make discretization appropriate for element, unless user specifically asks for fixed order
		int nDiscretizationPoint = (virtualRefinementOrder_ > 0) ? virtualRefinementOrder_ : curvOrder + 4;

		// Write element to VTK
		std::stringstream logstr;
		logstr << "CurvilinearVTKGridWriter: --Inserting entity " << gt;
		LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, logstr.str());
		baseWriter.template addCurvilinearElement<mydim, mydim>(gt, nodeSet, tagSet, curvOrder, nDiscretizationPoint, writeInterpolate_, writeExplode_);

		if (withField) {
			// Explicitly disallow writing fields for entities other than elements and faces
			assert((codim == ELEMENT_CODIM) || (codim == FACE_CODIM));
			writeField(baseWriter, fieldEntity);
		}
	}

	void writeField(BaseWriter & baseWriter, typename VTKVectorFunction3D::Entity fieldEntity)
	{
		writeScalarField(vtkScalarElementFunctionSet_, fieldEntity, baseWriter);
		writeVectorField(vtkVectorElementFunctionSet_, fieldEntity, baseWriter);
	}

	void writeField(BaseWriter & baseWriter, typename VTKVectorFunction2D::Entity fieldEntity)
	{
		writeScalarField(vtkScalarFaceFunctionSet_, fieldEntity, baseWriter);
		writeVectorField(vtkVectorFaceFunctionSet_, fieldEntity, baseWriter);
	}


	//template <class VTKVecPtrVector, class EntityType>
	template <class VTKFunction, class EntityType>
	void writeScalarField(std::vector<VTKFunction *> * vec, const EntityType & entity, BaseWriter & writer) {
		if (vec != nullptr) {
			for (unsigned int i = 0; i < vec->size(); i++)
			{
				LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Computing element scalar field" );
				(*vec)[i]->init(entity);

				LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Inserting element scalar field" );
				writer.template addScalarField<VTKFunction>(*((*vec)[i]));
			}
		}
	}


	template <class VTKFunction, class EntityType>
	void writeVectorField(std::vector<VTKFunction *> * vec, const EntityType & entity, BaseWriter & writer) {
		if (vec != nullptr) {
			for (unsigned int i = 0; i < vec->size(); i++)
			{
				LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Computing element vector field" );
				(*vec)[i]->init(entity);

				LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Inserting element vector field" );
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



private:

	const Grid & grid_;

	VTKScPtrVector2D  * vtkScalarFaceFunctionSet_;      // Store scalar functions defined over faces
	VTKVecPtrVector2D * vtkVectorFaceFunctionSet_;      // Store vector functions defined over faces
	VTKScPtrVector3D  * vtkScalarElementFunctionSet_;   // Store scalar functions defined over elements
	VTKVecPtrVector3D * vtkVectorElementFunctionSet_;   // Store vector functions defined over elements

	int virtualRefinementOrder_;   // User-defined discretization order for element sub-refinement
	bool writeDB_;       		// User can choose to write Domain Boundary Segments [Default - false]
	bool writeIB_;       		// User can choose to write Interior Boundary Segments [Default - false]
	bool writePB_;       		// User can choose to write Process Boundary Faces [Default - false]
	bool writePeriodic_;	// User can choose to write Periodic Boundary Faces [Default - false]
	bool writeGhost_;       		// User can choose to write Ghost Elements [Default - false]
	bool writeInterpolate_;       // User can choose to either interpolate a new grid over the element (true, recommended), or reuse the interpolatory points
	bool writeExplode_;            // User can choose to use explosion plotting for meshes
	bool writePatience_;           // User can choose to not write patience output for aestetical purposes

	bool writePeriodicBind_;   // User can choose to write additional info on the binding of neighboring periodic ghosts

	std::string vtkDataFormat_;  // Determines if VTU/PVTU files will use ascii or binary format for output. VTK is always ascii

};


} // namespace CurvGrid

} // Namespace Dune

#endif // DUNE_CURVILINEARVTKGRIDWRITER_HH
