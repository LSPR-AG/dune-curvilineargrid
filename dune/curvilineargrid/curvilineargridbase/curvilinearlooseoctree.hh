/***************************************************************************
                          CurvilinearLooseOctree.h  -  description
                             -------------------
    begin                : Wed Feb 25 2004
    copyright            : (C) 2004 by Roman Geus
    email                : roman.geus@psi.ch
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DUNE_CURVILINEARLOOSEOCTREE_HH
#define DUNE_CURVILINEARLOOSEOCTREE_HH

#include <vector>
#include <assert.h>

#include <dune/common/fvector.hh>


#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearoctant.hh>





namespace Dune {

/** 
 * Octree with overlapping octants
 * @param NodeType Type of nodes that are stored in the CurvilinearLooseOctree. Pointers to NodeType are stored
 *                 in the Octree to reference NodeType objects. NodeType must define a member function
 *                    NodeType::get_bounding_box(center, extent)
 *                 that computes the bounding box of a node and stores its center and its extent (the 
 *                 box sides halved) in center and extent respectively.
 * @author Roman Geus
 */
template <class ct, int cdim, class NodeType, class LogMsg>
class CurvilinearLooseOctree {
public:

	// typedefs
	// **************************************************************************
	typedef  ct      ctype;
	typedef  LogMsg  LoggingMessage;

	typedef Dune::FieldVector<ctype, cdim>                    Vertex;
	typedef Dune::CurvilinearOctant<ctype, cdim, NodeType>    CurvilinearOctant;

    // Logging Message Typedefs
    static const unsigned int LOG_CATEGORY_DEBUG = LoggingMessage::Category::DEBUG;


	/** Filter function deciding whether a point is inside a an OctreeNode */
	// FIXME: DO FUNCTOR DO NOT DO UGLY POINTER
	// FIXME: DO FUNCTOR DO NOT DO UGLY POINTER
	// FIXME: DO FUNCTOR DO NOT DO UGLY POINTER
	// FIXME: DO FUNCTOR DO NOT DO UGLY POINTER
	// FIXME: DO FUNCTOR DO NOT DO UGLY POINTER
	// FIXME: DO FUNCTOR DO NOT DO UGLY POINTER
	// FIXME: DO FUNCTOR DO NOT DO UGLY POINTER
	// FIXME: DO FUNCTOR DO NOT DO UGLY POINTER
	// FIXME: DO FUNCTOR DO NOT DO UGLY POINTER
	// FIXME: DO FUNCTOR DO NOT DO UGLY POINTER
    typedef bool (NodeType::*filter)(const Vertex&) const;

    typedef typename std::vector<NodeType *>::iterator NodePtrIterator;

public:

    /** Constructor: Initialise a CurvilinearLooseOctree living in a bounding box defined by
     * "center" and "length".
     */
    CurvilinearLooseOctree(
    		const Vertex& center,
    		double length,
    		int maxDepth,
    		MPIHelper & mpihelper) :
    			maxDepth_(maxDepth),
    			mpihelper_(mpihelper)
	{
    	root_ = new CurvilinearOctant(center, length);
	}


    /** Destructor: recursively deletes all Octants */
    ~CurvilinearLooseOctree()
    {
        delete root_;
        root_ = 0;
    }




    /** Add a new node, starting at the given Octree, and recursing at
     * max to level maxDepth_. The return value is a pointer to the
     * Octant storing the node.
     */
    CurvilinearOctant* addNode(NodeType* thisNode, CurvilinearOctant* octant=0, int depth=0)
    {
    	std::stringstream log_stream;
    	log_stream << "CurvilinearLooseOctree: Adding a node ElementIndex=" << thisNode->elementIndex() <<  " Octant=" << octant << " Depth=" << depth;

    	LoggingMessage::write<LOG_CATEGORY_DEBUG>(mpihelper_, __FILE__, __LINE__, log_stream.str());

        // root is the default octant
        if (octant == 0)  { octant = root_; }

        Vertex extent, center;
        thisNode->elementBoundingBox(center, extent);

        // if the octant is twice as big as the node,
        // we will add it to a child.
        if (
        	(depth < maxDepth_ ) &&
        	octant->length_ > 2*extent[0] &&
        	octant->length_ > 2*extent[1] &&
        	octant->length_ > 2*extent[2]
           )
        {
            // find child octant containing the mid-point of box
            int childIndex = octant->childIndex(center);

            // If Octant does not yet exist: create new one
            if (octant->children_[childIndex] == 0)
            	octant->children_[childIndex] = new CurvilinearOctant(octant, childIndex);

            // call myself recursively
            return addNode(thisNode, octant->children_[childIndex], depth + 1);
        } else
        {
            // found destination octant: add node
        	octant->addNode(thisNode);
            return octant;
        }
    }


    /** Find nodes who contain the given point. For determining
        whether a node contains the point the "filter" function is
        called. All found nodes are appended to
        "node_list". "nNodeVisited" is a variable which is incremented
        each time a node is checked. On input "nNodeVisited" should be
        set to zero. "octant" specifies the starting octant. By
        default the search starts at the root node.
    */
    void findNode(const Vertex& coord,
                             std::vector<int>& elementIndices,
                             int& nNodeVisited,
                             CurvilinearLooseOctree::filter filter,
                             CurvilinearOctant* octant=0)
    {
        // root is the default octant
        if (octant == 0)  { octant = root_; }

        // return if coord is outside octant (including overlapping region)
        Vertex dist = coord - octant->_center;
        double extentLength = 2.0 * octant->length_;
        if (fabs(dist[0]) > extentLength ||
            fabs(dist[1]) > extentLength ||
            fabs(dist[2]) > extentLength)
        { return; }

        // check nodes stored in current Octant
        NodePtrIterator it = octant->node_.begin();
        while (it != octant->node_.end()) {
            NodeType* thisNode = *it;
            if ((thisNode->*filter)(coord))  { elementIndices.push_back(thisNode->elementIndex()); }
            nNodeVisited++;
            it++;
        }

        // descend to children
        for (int i = 0; i < 8; i ++)
        {
            if (octant->children_[i])  { findNode(coord, elementIndices, nNodeVisited, filter, octant->children_[i]); }
        }


    }


    /** Finds a node who contain the given point. For determining
        whether a node contains the point the "filter" function is
        called. "nNodeVisited" is a variable which is incremented each
        time a node is checked. On input "nNodeVisited" should be set
        to zero. "octant" specifies the starting octant. By default
        the search starts at the root node. If no node was found, 0 is
        returned.
    */
    int findSingleNode(const Vertex& coord,
                                       int& nNodeVisited,
                                       CurvilinearLooseOctree::filter filter,
                                       CurvilinearOctant* octant=0)
    {
        NodeType* thisNode;

        // root is the default octant
        if (octant == 0)  { octant = root_; }

        // return if coord is outside octant (including overlapping region)
        Vertex dist = coord - octant->_center;
        double extentLength = 2.0 * octant->length_;
        if (fabs(dist[0]) > extentLength ||
            fabs(dist[1]) > extentLength ||
            fabs(dist[2]) > extentLength)
        { return 0; }

        // check nodes stored in current Octant
        NodePtrIterator it = octant->node_.begin();
        while (it != octant->node_.end()) {
            thisNode = *it;
            nNodeVisited ++;
            it ++;
            if ((thisNode->*filter)(coord))  { return thisNode->elementIndex(); }
        }
        // descend to children
        for (int i = 0; i < 8; i ++) {
            if (octant->children_[i]) {
                thisNode = findSingleNode(coord, nNodeVisited, filter, octant->children_[i]);

                // if a node was found abort descend to other children
                if (thisNode)  { return thisNode; }
            }
        }
        // Nothing found, return 0
        return 0;
    }


    /** Compute statistics on the octree (maximum depth over all octants, average
        depth over all nodes, number of octants, number of nodes).
    */
    void statistics(int& maxDepth,
                        double& avgNodeDepth,
                        int& nOctant,
                        int& nNode)
    {
        CurvilinearOctant* octant = root_;
        int depth = 0;
        double sumDepth = 0.0;

        maxDepth = 0;
        nOctant = 0;
        nNode = 0;
        statisticsRecursive(octant, depth, maxDepth, sumDepth, nOctant, nNode);
        avgNodeDepth = sumDepth / nNode;
    }



    /** Export CurvilinearLooseOctree to a series of VTK files (one for each level).
        The VTK files are named "octree_lX.vtk", where X is the level
        number. The root Octant has level 0.
     */
    void vtkWrite()
    {
        // get max depth
        int maxDepth, nOctant, nNode;
        double avgNodeDepth;

        statistics(maxDepth, avgNodeDepth, nOctant, nNode);

        // export all levels to separate VTK files
        for (int level = 0; level <= maxDepth; ++level) {
            std::ostringstream buf;
            buf << "octree_l" << level << ".vtk";
            vtkWriteLevel(level, buf.str());
        }
    }


    /** Export one level of the CurvilinearLooseOctree to a VTK file.
     */
    void vtkWriteLevel(int level, std::string filename)
    {
        std::vector<int> cells;
        Vertex center;
        center[0] = 0;
        center[1] = 1;
        center[2] = 2;

        // recursive descend
        vtkDataRecursive(root_, 0, center, level, cells);
        int nOctant = cells.size() / 9;

        std::ofstream of;
        of.open(filename.c_str());
        //rAssert(of.is_open());

        // VTK header
        of.precision(6);
        //of << scientific;
        of << "# vtk DataFile Version 2.0" << std::endl;
        of << "generated using VtkExport::export_eigenfields" << std::endl;
        of << "ASCII" << std::endl << std::endl;
        of << "DATASET UNSTRUCTURED_GRID" << std::endl;

        // export all points in level, (2**depth + 1)**3 points
        int n = (1 << level) + 1;
        Vertex origin = root_->center_;
        origin -= root_->length_;
        double side = 2.0*root_->length_ / (1 << level);
        of << "POINTS " << n*n*n << " float" << std::endl;

        for (int iterZ = 0; iterZ < n; ++ iterZ)
        {
            for (int iterY = 0; iterY < n; ++ iterY)
            {
                for (int iterX = 0; iterX < n; ++ iterX)
                {
                    of << origin[0] + iterX*side << " "
                       << origin[1] + iterY*side << " "
                       << origin[2] + iterZ*side << std::endl;
                }
            }
        }
        of << std::endl;

        // cells
        of << "CELLS " << nOctant << " " << cells.size() << std::endl;
        std::vector<int>::iterator it = cells.begin();
        for (int i = 0; it != cells.end(); ++ it, ++ i) {
            of << *it;
            if (i % 9 == 8)  { of << std::endl; }
            else             { of << " "; }
        }
        of << std::endl;

        // cell types
        of << "CELL_TYPES " << nOctant << std::endl;
        for (int i = 0; i < nOctant; ++ i)  { of << "11" << std::endl; }
        of << std::endl;

        // cell data
        of << "CELL_DATA " << nOctant << std::endl;
        of << "SCALARS cell_scalars int 1" << std::endl;
        of << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < nOctant; ++ i)  { of << "0" << std::endl; }
        of << std::endl;
    }


protected:

    void statisticsRecursive(CurvilinearOctant* octant,
                                  int depth,
                                  int& maxDepth,
                                  double& sumDepth,
                                  int& nOctant,
                                  int& nNode)
    {
        // update statistics
        if (depth > maxDepth)  { maxDepth = depth; }
        sumDepth += octant->node_.size() * depth;
        nOctant ++;
        nNode += octant->node_.size();

        // descend to children
        for (int i = 0; i < 8; i ++)
        {
            if (octant->children_[i])  { statisticsRecursive(octant->children_[i], depth+1, maxDepth, sumDepth, nOctant, nNode); }
        }
    }


    void vtkDataRecursive(CurvilinearOctant* octant,
    		                  int depth,
                              Vertex & center,
                              int targetDepth,
                              std::vector<int>& cells)
    {
        if (depth == targetDepth)
        {
            // export octant
            int n = (1 << depth) + 1;
            cells.push_back(8);
            cells.push_back((center[0]+0) + (center[1]+0)*n + (center[2]+0)*n*n);
            cells.push_back((center[0]+1) + (center[1]+0)*n + (center[2]+0)*n*n);
            cells.push_back((center[0]+0) + (center[1]+1)*n + (center[2]+0)*n*n);
            cells.push_back((center[0]+1) + (center[1]+1)*n + (center[2]+0)*n*n);
            cells.push_back((center[0]+0) + (center[1]+0)*n + (center[2]+1)*n*n);
            cells.push_back((center[0]+1) + (center[1]+0)*n + (center[2]+1)*n*n);
            cells.push_back((center[0]+0) + (center[1]+1)*n + (center[2]+1)*n*n);
            cells.push_back((center[0]+1) + (center[1]+1)*n + (center[2]+1)*n*n);
        } else {
            // descend to children
            for (int i = 0; i < 8; i ++) {
                if (octant->children_[i])
                {
                	Vertex subCenter;
                	subCenter[0] = center[0] + ((i & 0x1) == 0x1) * (1 << (targetDepth - depth - 1));
                	subCenter[1] = center[1] + ((i & 0x2) == 0x2) * (1 << (targetDepth - depth - 1));
                	subCenter[2] = center[2] + ((i & 0x4) == 0x4) * (1 << (targetDepth - depth - 1));

                    vtkDataRecursive(octant->children_[i], depth+1, subCenter, targetDepth, cells);
                }
            }
        }
    }

private:

    MPIHelper &mpihelper_;

    CurvilinearOctant* root_;    /** Root Octant */
    int maxDepth_;               /** Maximum tree depth */
};


} // namespace Dune

#endif // DUNE_CURVILINEARLOOSEOCTREE_HH
