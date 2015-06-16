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
	typedef typename GridType::GridBaseType::template Codim<ELEMENT_CODIM>::EntityGeometry   EntityGeometry;
	typedef typename EntityGeometry::GlobalCoordinate  GlobalCoordinate;

	typedef typename Grid::PhysicalTagType           PhysicalTagType;
	typedef typename Grid::InterpolatoryOrderType    InterpolatoryOrderType;

	typedef typename GridType::LoggingMessage           LoggingMessage;
	static const unsigned int LOG_CATEGORY_DEBUG = LoggingMessage::Category::DEBUG;


	typedef Dune::VTKFunction<LeafGridView>   VTKFunction;




public:
	CurvilinearVTKGridWriter(const Grid & grid) :
		grid_(grid),
		writer_(
			grid.comm().rank(),
			grid.comm().size()
		)
	{

	}


	// Write Grid only to VTK
	void write(std::string filename)
	{
		std::vector<VTKFunction *> emptyFieldSet;
		write(filename, emptyFieldSet);
	}


	// Write Grid to VTK, including vector field set defined over each element
	void write(std::string filename, std::vector<VTKFunction *> & vtkFunctionSet)
	{
		// Properties
        int nDiscretizationPoint = 6;
        bool interpolate = true;
        bool explode = false;
        std::vector<bool>  writeCodim { true, true, false, false };  // Use elements for discretization, and do not use entities of other codimensions

		// Iterate over elements
  		/** \brief Iterate ove all elements of Interior Border partition */
		LeafGridView leafView = grid_.leafGridView();

		LoggingMessage::template writeStatic<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Writing Elements" );
		for (auto&& elementThis : elements(leafView, Dune::Partitions::all))
		{
			Dune::GeometryType geomtype            = elementThis.type();
			Dune::PartitionType thisPType          = elementThis.partitionType();

			// Constructing a geometry is quite expensive, do it only once
			//EntityGeometry geom = it->geometry();
			EntityGeometry geom                    = grid_.template entityBaseGeometry<ELEMENT_CODIM>(elementThis);
			std::vector<GlobalCoordinate>  nodeSet = geom.vertexSet();

			PhysicalTagType         physicalTag  = grid_.template entityPhysicalTag<ELEMENT_CODIM>(elementThis);
			InterpolatoryOrderType  elementOrder = grid_.template entityInterpolationOrder<ELEMENT_CODIM>(elementThis);
			std::vector<int>        tagSet  { physicalTag, thisPType, grid_.comm().rank() };

			// Write element to VTK
			LoggingMessage::template writeStatic<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: --Inserting element geometry" );
			writer_.template addCurvilinearElement<dimension - ELEMENT_CODIM>(geomtype, nodeSet, tagSet, elementOrder, nDiscretizationPoint, interpolate, explode, writeCodim);


			if (thisPType != Dune::PartitionType::GhostEntity)
			{
				// Write faces - Domain and Process Boundaries
				for (auto&& intersection : intersections(leafView, elementThis))
				{
					const EntityFace faceThis = elementThis.template subEntity<FACE_CODIM>(intersection.indexInInside());

					Dune::GeometryType  faceGeomType   = faceThis.type();
					Dune::PartitionType thisFacePType  = faceThis.partitionType();

					if ((thisFacePType == Dune::PartitionType::BorderEntity) || intersection.boundary())
					{
						// Constructing a geometry is quite expensive, do it only once
						//EntityGeometry geom = it->geometry();
						EntityGeometry facegeom                    = grid_.template entityBaseGeometry<ELEMENT_CODIM>(faceThis);
						std::vector<GlobalCoordinate>  faceNodeSet = facegeom.vertexSet();

						PhysicalTagType         faceTag            = grid_.template entityPhysicalTag<ELEMENT_CODIM>(faceThis);
						std::vector<int>        faceTagSet  { faceTag, thisFacePType, grid_.comm().rank() };

						// Write element to VTK
						LoggingMessage::template writeStatic<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: --Inserting face geometry" );
						writer_.template addCurvilinearElement<dimension - FACE_CODIM>(faceGeomType, faceNodeSet, faceTagSet, elementOrder, nDiscretizationPoint, interpolate, explode, writeCodim);
					}
				}


				// It is unexpected that the process would know the field on its own Ghost element. It is most likely that for memory reasons this field
				// will be stored only once, and on a process for which this element is an interior element
				for (int i = 0; i < vtkFunctionSet.size(); i++)
				{
					LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Computing element field" );
					// Construct VTKElementaryFunction
					VTKElementaryFunction<EntityElement, VTKFunction> vtkefunc(elementThis, vtkFunctionSet[i]);

					LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: ---Inserting element field" );
					// Write field
					writer_.addField(vtkFunctionSet[i]->name(), vtkefunc);
				}
			}
		}


		// Write all faces
		LoggingMessage::template writeStatic<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Writing Faces" );
		for (auto&& faceThis : faces(leafView, Dune::Partitions::interiorBorder))
		{

		}


		// Write the data to vtk file
		LoggingMessage::template writeStatic<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Writing .vtk file = " + filename);
		writer_.writeVTK(filename);

		// Delete vtk functions
		for (int i = 0; i < vtkFunctionSet.size(); i++)  { delete vtkFunctionSet[i]; }
	}








private:

	const Grid & grid_;
	BaseWriter writer_;

};




} // Namespace Dune

#endif // DUNE_CURVILINEARVTKGRIDWRITER_HH
