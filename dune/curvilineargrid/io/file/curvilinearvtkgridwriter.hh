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

	typedef typename GridType::template Codim<ELEMENT_CODIM>::Entity                         EntityElement;
	typedef typename GridType::template Codim<FACE_CODIM>::Entity                            EntityFace;
	typedef typename GridType::LeafGridView::Traits::Intersection  Intersection;
	typedef typename GridType::GridBaseType::template Codim<ELEMENT_CODIM>::EntityGeometry   ElementGeometry;
	typedef typename GridType::GridBaseType::template Codim<FACE_CODIM>::EntityGeometry      FaceGeometry;
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
		baseWriter_(grid_.comm().rank(), grid_.comm().size()),
		vtkScalarFaceFunctionSet_(nullptr),
		vtkScalarElementFunctionSet_(nullptr),
		vtkVectorFaceFunctionSet_(nullptr),
		vtkVectorElementFunctionSet_(nullptr),
		virtualRefinementOrder_(0),
		writePB_(false),
		writeDB_(false),
		writeIB_(false),
		writePeriodic_(false),
		writeGhost_(false),
		writeInterpolate_(true),
		writeExplode_(false),
		writePatience_(true)
	{
		// Use elements and faces for discretization, and do not use entities of other codimensions,
		// as it is an overkill for visualisation of solutions over a functional mesh
        writeCodim_ = std::vector<bool> { true, true, false, false };
	}


	// Allow user to fix virtual refiniment level
	void useFixedVirtualRefinement(int newOrder)  { virtualRefinementOrder_ = newOrder; }

	// Allow user to only write interior elements to accelerate writer
	void writeProcessBoundary(bool write)  { writePB_ = write; }
	void writeDomainBoundary(bool write)  { writeDB_ = write; }
	void writeInteriorBoundary(bool write)  { writeIB_ = write; }
	void writePeriodicBoundary(bool write)  { writePeriodic_ = write; }
	void writeGhost(bool write)						{ writeGhost_ = write; }

	// Allow user to only write interior elements to accelerate writer
	void writePatience(bool patience)  { writePatience_ = patience; }

	// User can choose to use explosion plotting for meshes
	void writeInterpolate(bool interpolate)  { writeInterpolate_ = interpolate; }
	void writeExplode(bool explode)  { writeExplode_ = explode; }

	// Add a scalar or a vector field set to the grid
	// NOTE: IT IS VERY IMPORTANT THAT ALL FIELDS HAVE A UNIQUE NAME
	void addFieldSet(VTKScPtrVector2D & field)   { vtkScalarFaceFunctionSet_    = &field; }
	void addFieldSet(VTKScPtrVector3D & field)   { vtkScalarElementFunctionSet_ = &field; }
	void addFieldSet(VTKVecPtrVector2D & field)  { vtkVectorFaceFunctionSet_    = &field; }
	void addFieldSet(VTKVecPtrVector3D & field)  { vtkVectorElementFunctionSet_ = &field; }


	// Write Grid to VTK, including vector field set defined over each element
	void write(std::string path, std::string filenamePrefix)
	{
		// Determine Grid Capabilities
		bool isGridParallel = grid_.comm().size() > 0;
		bool hasGridGhost = grid_.ghostSize(ELEMENT_CODIM) > 0;

        // Determine if any of the faces will be written
        bool writeAnyFace = writeDB_ || writePB_ || writeIB_ || writePeriodic_;

        // Initialize all vector fields, in case there are none on this process
		if (vtkScalarFaceFunctionSet_ != nullptr)  { for (unsigned int i = 0; i < vtkScalarFaceFunctionSet_->size(); i++) { baseWriter_.initScalarField( *((*vtkScalarFaceFunctionSet_)[i])); } }
		if (vtkVectorFaceFunctionSet_ != nullptr)  { for (unsigned int i = 0; i < vtkVectorFaceFunctionSet_->size(); i++) { baseWriter_.initVectorField( *((*vtkVectorFaceFunctionSet_)[i])); } }
		if (vtkScalarElementFunctionSet_ != nullptr)  { for (unsigned int i = 0; i < vtkScalarElementFunctionSet_->size(); i++) { baseWriter_.initScalarField( *((*vtkScalarElementFunctionSet_)[i])); } }
		if (vtkVectorElementFunctionSet_ != nullptr)  { for (unsigned int i = 0; i < vtkVectorElementFunctionSet_->size(); i++) { baseWriter_.initVectorField( *((*vtkVectorElementFunctionSet_)[i])); } }

		// Iterate over elements
  		/** \brief Iterate over all elements of Interior Border partition */
        int elemIterCount = 0;         // Count elements
		LeafGridView leafView = grid_.leafGridView();
		LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Writing Elements" );
		for (auto&& elementThis : elements(leafView, Dune::Partitions::interiorBorder))
		{
			if (writePatience_) {
				LoggingMessage::writePatience("CurvilinearVTKGridWriter: Writing Elements to VTK...", elemIterCount++, grid_.size(ELEMENT_CODIM));
			}

			// Write this interior element
			Dune::PartitionType elementType = elementThis.partitionType();
			bool interiorElementWriteField = true;
			writeEntity<ELEMENT_CODIM, EntityElement, EntityElement>(elementThis, elementThis, elementType, interiorElementWriteField);

			// Write faces and ghost elements
			if (writeAnyFace || writeGhost_)
			{
				for (auto&& intersection : intersections(leafView, elementThis))
				{
					// Get geometric parameters
					const EntityFace faceThis = elementThis.template subEntity<FACE_CODIM>(intersection.indexInInside());
					PhysicalTagType faceTag = grid_.template entityPhysicalTag<FACE_CODIM>(faceThis);
					Dune::GeometryType  faceGeomType   = faceThis.type();

					// NOTE: CurvGridVTKWriter introduces extended partition type format to better visualise special boundary segments
					int thisFaceEffectivePartitionType  = faceThis.partitionType();

					// Determine if the face will be written, based on its type
					// Note that periodic boundaries are also domain boundaries
					bool writeFace = false;
					bool isPeriodic = ((intersection.neighbor() == true) && (intersection.boundary() == true));
					bool isDB = ((intersection.neighbor() == false) && (intersection.boundary() == true)) || isPeriodic;
					bool isPB = isGridParallel
							? ((intersection.neighbor() == true) && (intersection.outside().partitionType() == Dune::PartitionType::GhostEntity))
							: ((intersection.neighbor() == false) && (intersection.boundary() == false));
					bool isIB = ((!isDB) && (faceTag >= 0));

					if				( isPeriodic && writePeriodic_)	{ writeFace = true;  thisFaceEffectivePartitionType = PERIODIC_BOUNDARY_PARTITION_TYPE; }
					else if		( isDB && writeDB_)					{ writeFace = true;  thisFaceEffectivePartitionType = BOUNDARY_SEGMENT_PARTITION_TYPE; }
					else if		( isIB && writeIB_)						{ writeFace = true;  thisFaceEffectivePartitionType = INTERIOR_BOUNDARY_SEGMENT_PARTITION_TYPE; }
					else if		( isPB && writePB_)						{
						// Do not write PB twice, it is wasteful. Determine which process writes it using unit normal at center
						// [TODO] Normal orientation is ugly. Can determine owner natively, if can extract easily neighbor rank
						GlobalCoordinate faceNormal = intersection.unitOuterNormal( ReferenceElements2d::simplex().position( 0, 0 ) );
						writeFace = (faceNormal[0] > 0) ||
								((faceNormal[0] == 0) && (faceNormal[1] > 0)) ||
								((faceNormal[0] == 0) && (faceNormal[1] == 0) && (faceNormal[2] > 0));
					}

					// Write face
					// Note: For now, only write face-based fields for domain boundaries
					if  (writeFace) { writeEntity<FACE_CODIM, EntityFace, Intersection>(faceThis, intersection, thisFaceEffectivePartitionType, isDB); }

					// Write ghost element
					if (writeGhost_ && (isPB || isPeriodic)) {
						EntityElement thisGhost = intersection.outside();
						Dune::PartitionType ghostType = thisGhost.partitionType();

						if (isPB) { assert(ghostType == Dune::PartitionType::GhostEntity); }
						int ghostEffectivePartitionType = isPB ? ghostType : PERIODIC_GHOST_PARTITION_TYPE;

						// Note: It is unexpected that the process would know the field on its own Ghost element.
						// It is most likely that for memory reasons that this field will be stored only once, and on a process for which this element is an interior element
						bool ghostElementWriteField = true;
						writeEntity<ELEMENT_CODIM, EntityElement, EntityElement>(thisGhost, thisGhost, ghostEffectivePartitionType, ghostElementWriteField);
					}
				}
			}
		}

		// Write the data to vtk file
		LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Writing .vtk file = " + path + filenamePrefix + ".pvtu");

		//writer_.writeVTK(filenamePrefix + ".vtk");
		baseWriter_.writeParallelVTU(path, filenamePrefix);

		// With current design it is expected that the VTK Writer deletes all the VTK functions
		deleteField(vtkScalarFaceFunctionSet_);
		deleteField(vtkVectorFaceFunctionSet_);
		deleteField(vtkScalarElementFunctionSet_);
		deleteField(vtkVectorElementFunctionSet_);
	}




protected:

	template <int codim, class GridEntityType, class FieldEntityType>
	void writeEntity(const GridEntityType & gridEntity, const FieldEntityType & fieldEntity, int ptEff, bool withField) {
		const int mydim = dimension - codim;
		Dune::GeometryType gt   = gridEntity.type();

		// Constructing a geometry is quite expensive, do it only once
		typedef typename GridType::GridBaseType::template Codim<codim>::EntityGeometry   ThisGeometry;
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
		baseWriter_.template addCurvilinearElement<mydim>(gt, nodeSet, tagSet, curvOrder, nDiscretizationPoint, writeInterpolate_, writeExplode_, writeCodim_);

		if (withField) {
			// Explicitly disallow writing fields for entities other than elements and faces
			assert((codim == ELEMENT_CODIM) || (codim == FACE_CODIM));
			writeField(fieldEntity);
		}
	}

	void writeField(typename VTKVectorFunction<Grid, DIM3D-ELEMENT_CODIM>::Entity fieldEntity)
	{
		writeScalarField(vtkScalarElementFunctionSet_, fieldEntity, baseWriter_);
		writeVectorField(vtkVectorElementFunctionSet_, fieldEntity, baseWriter_);
	}

	void writeField(typename VTKVectorFunction<Grid, DIM3D-FACE_CODIM>::Entity fieldEntity)
	{
		writeScalarField(vtkScalarFaceFunctionSet_, fieldEntity, baseWriter_);
		writeVectorField(vtkVectorFaceFunctionSet_, fieldEntity, baseWriter_);
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
	BaseWriter baseWriter_;

	VTKScPtrVector2D  * vtkScalarFaceFunctionSet_;      // Store scalar functions defined over faces
	VTKVecPtrVector2D * vtkVectorFaceFunctionSet_;      // Store vector functions defined over faces
	VTKScPtrVector3D  * vtkScalarElementFunctionSet_;   // Store scalar functions defined over elements
	VTKVecPtrVector3D * vtkVectorElementFunctionSet_;   // Store vector functions defined over elements

	int virtualRefinementOrder_;   // User-defined discretization order for element sub-refinement
	std::vector<bool> writeCodim_; // Determine the codimensions of entities to be written
	bool writeDB_;       		// User can choose to write Domain Boundary Segments [Default - false]
	bool writeIB_;       		// User can choose to write Interior Boundary Segments [Default - false]
	bool writePB_;       		// User can choose to write Process Boundary Faces [Default - false]
	bool writePeriodic_;	// User can choose to write Periodic Boundary Faces [Default - false]
	bool writeGhost_;       		// User can choose to write Ghost Elements [Default - false]
	bool writeInterpolate_;       // User can choose to either interpolate a new grid over the element (true, recommended), or reuse the interpolatory points
	bool writeExplode_;            // User can choose to use explosion plotting for meshes
	bool writePatience_;           // User can choose to not write patience output for aestetical purposes

};


} // namespace CurvGrid

} // Namespace Dune

#endif // DUNE_CURVILINEARVTKGRIDWRITER_HH
