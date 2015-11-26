#ifndef DUNE_CURVILINEARVTKGRIDWRITER_HH
#define DUNE_CURVILINEARVTKGRIDWRITER_HH

#include <dune/curvilineargrid/io/file/curvilinearvtkwriter.hh>

/** \brief This is only needed to define the VTK Function at the moment **/
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

namespace Dune
{


template <class Element, class VTKFunction>
class VTKElementaryFunction
{
	static const int dimension = Element::dimension;

	typedef Dune::FieldVector<double, dimension>  LocalCoordinate;
	typedef Dune::FieldVector<double, dimension>  GlobalCoordinate;

public:

	VTKElementaryFunction(const Element & element, const VTKFunction * vtkfunction) :
		element_(element),
		vtkfunction_(vtkfunction)
	{

	}


	GlobalCoordinate evaluate(const LocalCoordinate &xi)
	{
		GlobalCoordinate rez;
		for (int i = 0; i < dimension; i++)  { rez[i] = vtkfunction_->evaluate(i, element_, xi); }
		return rez;
	}


private:
	const Element & element_;
	const VTKFunction * vtkfunction_;
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

	typedef Dune::VTKFunction<LeafGridView>   VTKFunction;




public:
	CurvilinearVTKGridWriter(const Grid & grid) :
		grid_(grid),
		virtualRefinementOrder_(0),
		writeInteriorOnly_(false),
		writePatience_(true)
	{

	}


	// Allow user to fix virtual refiniment level
	void useFixedVirtualRefinement(int newOrder)  { virtualRefinementOrder_ = newOrder; }

	// Allow user to only write interior elements to accelerate writer
	void writeInteriorOnly(bool interiorOnly)  { writeInteriorOnly_ = interiorOnly; }

	// Allow user to only write interior elements to accelerate writer
	void writePatience(bool patience)  { writePatience_ = patience; }

	// Write Grid only to VTK
	void write(std::string path, std::string filenamePrefix)
	{
		std::vector<VTKFunction *> emptyFieldSet;
		write(path, filenamePrefix, emptyFieldSet);
	}


	// Write Grid to VTK, including vector field set defined over each element
	void write(std::string path, std::string filenamePrefix, std::vector<VTKFunction *> & vtkFunctionSet)
	{
		// Declare the curvilinear vtk writer
		BaseWriter writer(grid_.comm().rank(), grid_.comm().size());

		// Properties
        bool interpolate = true;
        bool explode = false;
        std::vector<bool>  writeCodim { true, true, false, false };  // Use elements for discretization, and do not use entities of other codimensions

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

			if ((!writeInteriorOnly_) || (thisPType != Dune::PartitionType::GhostEntity))
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
				writer.template addCurvilinearElement<dimension - ELEMENT_CODIM>(geomtype, nodeSet, tagSet, elementOrder, nDiscretizationPoint, interpolate, explode, writeCodim);


				if (thisPType != Dune::PartitionType::GhostEntity)
				{
					// It is unexpected that the process would know the field on its own Ghost element. It is most likely that for memory reasons this field
					// will be stored only once, and on a process for which this element is an interior element
					for (unsigned int i = 0; i < vtkFunctionSet.size(); i++)
					{

						LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Computing element field" );
						// Construct VTKElementaryFunction
						VTKElementaryFunction<EntityElement, VTKFunction> vtkefunc(elementThis, vtkFunctionSet[i]);

						LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Inserting element field" );
						// Write field
						writer.addField(vtkFunctionSet[i]->name(), vtkefunc);
					}

					// Write faces - Domain and Process Boundaries
					if (!writeInteriorOnly_) // Do it only if user wants to see them
					{
						for (auto&& intersection : intersections(leafView, elementThis))
						{
							const EntityFace faceThis = elementThis.template subEntity<FACE_CODIM>(intersection.indexInInside());

							Dune::GeometryType  faceGeomType   = faceThis.type();
							StructuralType      thisFacePType  = faceThis.partitionType();

							if ((thisFacePType == Dune::PartitionType::BorderEntity) || intersection.boundary())
							{
								if (intersection.boundary())  { thisFacePType = CurvGrid::BOUNDARY_SEGMENT_PARTITION_TYPE; }

								// Constructing a geometry is quite expensive, do it only once
								//EntityGeometry geom = it->geometry();
								FaceGeometry facegeom                    = grid_.template entityBaseGeometry<FACE_CODIM>(faceThis);
								std::vector<GlobalCoordinate>  faceNodeSet = facegeom.vertexSet();

								PhysicalTagType         faceTag          = grid_.template entityPhysicalTag<FACE_CODIM>(faceThis);
								std::vector<int>        faceTagSet  { faceTag, thisFacePType, grid_.comm().rank() };

								// Write element to VTK
								LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: --Inserting face geometry" );
								writer.template addCurvilinearElement<dimension - FACE_CODIM>(faceGeomType, faceNodeSet, faceTagSet, elementOrder, nDiscretizationPoint, interpolate, explode, writeCodim);
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

		// Delete vtk functions
		for (unsigned int i = 0; i < vtkFunctionSet.size(); i++)  { delete vtkFunctionSet[i]; }
		vtkFunctionSet.clear();

	}


private:

	const Grid & grid_;
	int virtualRefinementOrder_;   // User-defined discretization order for element sub-refinement
	bool writeInteriorOnly_;       // User can choose to only write interior elements, thus accelerating writing procedure
	bool writePatience_;           // User can choose to not write patience output for aestetical purposes

};




} // Namespace Dune

#endif // DUNE_CURVILINEARVTKGRIDWRITER_HH
