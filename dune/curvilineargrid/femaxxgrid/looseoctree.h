/***************************************************************************
                          looseoctree.h  -  description
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

#ifndef LOOSEOCTREE_H
#define LOOSEOCTREE_H

#include <vector>
#include "octant.h"

using namespace std;

namespace mesh {

/** 
 * Octree with overlapping octants
 * @param NodeType Type of nodes that are stored in the LooseOctree. Pointers to NodeType are stored 
 *                 in the Octree to reference NodeType objects. NodeType must define a member function
 *                    NodeType::get_bounding_box(node_center, node_extent)
 *                 that computes the bounding box of a node and stores its center and its extent (the 
 *                 box sides halved) in node_center and node_extent respectively.
 * @author Roman Geus
 */
template <typename NodeType>
class LooseOctree {
public: // typedefs
    /** Filter function deciding whether a point is inside a an OctreeNode
     */
    typedef bool (NodeType::*filter)(const Vector3&) const;
public:
    /** Constructor: Initialise a LooseOctree living in a bounding box defined by
     * "center" and "length".
     */
    LooseOctree(const Vector3& center, double length, int max_depth);
    /** Destructor: recursively deletes all Octants
     */
    ~LooseOctree();
    /** Add a new node, starting at the given Octree, and recursing at
     * max to level _max_depth. The return value is a pointer to the
     * Octant storing the node.
     */
    Octant<NodeType>* add_node(NodeType* node, Octant<NodeType>* start=0, int depth=0);
    /** Find nodes who contain the given point. For determining
        whether a node contains the point the "filter" function is
        called. All found nodes are appended to
        "node_list". "nof_visited" is a variable which is incremented
        each time a node is checked. On input "nof_visited" should be
        set to zero. "octant" specifies the starting octant. By
        default the search starts at the root node.
    */
    void find_nodes_by_point(const Vector3& coord,
                             std::vector<NodeType *>& node_list,
                             int& nof_visited,
                             LooseOctree::filter filter,
                             Octant<NodeType>* octant=0);
    /** Finds a node who contain the given point. For determining
        whether a node contains the point the "filter" function is
        called. "nof_visited" is a variable which is incremented each
        time a node is checked. On input "nof_visited" should be set
        to zero. "octant" specifies the starting octant. By default
        the search starts at the root node. If no node was found, 0 is
        returned.
    */
    NodeType* find_one_node_by_point(const Vector3& coord,
                                       int& nof_visited,
                                       LooseOctree::filter filter,
                                       Octant<NodeType>* octant=0);
     
    /** Compute statistics on the octree (maximum depth over all octants, average
        depth over all nodes, number of octants, number of nodes).
    */
    void get_statistics(int& max_depth,
                        double& avg_node_depth,
                        int& nof_octants,
                        int& nof_nodes);

    /** Export LooseOctree to a series of VTK files (one for each level).

        The VTK files are named "octree_lX.vtk", where X is the level
        number. The root Octant has level 0.
     */
    void export_vtk();
    /** Export one level of the LooseOctree to a VTK file.
     */
    void export_vtk_level(int level, std::string filename);

protected:
    void get_statistics_recursive(Octant<NodeType>* octant,
                                  int depth,
                                  int& max_depth,
                                  double& sum_depth,
                                  int& nof_octants,
                                  int& nof_nodes);

    void export_vtk_recursive(Octant<NodeType>* octant,int depth,
                              int x,
                              int y,
                              int z,
                              int target_depth,
                              std::vector<int>& cells);

private:
    /** Root Octant */
    Octant<NodeType>* _root;
    /** Maximum tree depth */
    int _max_depth;
};
  

template <typename NodeType>
LooseOctree<NodeType>::LooseOctree(const Vector3& center, double length, int max_depth)
    : _max_depth(max_depth)
{
    _root = new Octant<NodeType>(center, length);
}

template <typename NodeType>
LooseOctree<NodeType>::~LooseOctree() {
    delete _root;
    _root = 0;
}

template <typename NodeType>
Octant<NodeType>* LooseOctree<NodeType>::add_node(NodeType* node, Octant<NodeType>* octant, int depth) {
    // root is the default octant
    if (octant == 0)
        octant = _root;
      
    Vector3 node_extent, node_center;
    node->get_bounding_box(node_center, node_extent);

    // if the octant is twice as big as the node,
    // we will add it to a child.
    if ((depth < _max_depth ) &&
        octant->_length > 2*node_extent.x &&
        octant->_length > 2*node_extent.y &&
        octant->_length > 2*node_extent.z
        ) {
        // find child octant containing the mid-point of box
        int child_idx = octant->get_child_index_for_point(node_center);

        // If Octant does not yet exist: create new one
        if (octant->_children[child_idx] == 0)
            octant->_children[child_idx] = new Octant<NodeType>(octant, child_idx);

        // call myself recursively
        return add_node(node, octant->_children[child_idx], depth + 1);

    } else {
        // found destination octant: add node
        octant->add_node(node);
        return octant;
    }
}

template <typename NodeType>
void LooseOctree<NodeType>::find_nodes_by_point(const Vector3& coord,
                                      std::vector<NodeType *>& nodes,
                                      int& nof_visited,
                                      LooseOctree<NodeType>::filter filter,
                                      Octant<NodeType>* octant) {
    // root is the default octant
    if (octant == 0)
        octant = _root;
    // return if coord is outside octant (including overlapping region)
    Vector3 dist = coord - octant->_center;
    double ext_length = 2.0 * octant->_length;
    if (fabs(dist.x) > ext_length ||
        fabs(dist.y) > ext_length ||
        fabs(dist.z) > ext_length)
        return;
    // check nodes stored in current Octant
    typename std::vector<NodeType *>::iterator it = octant->_nodes.begin();
    while (it != octant->_nodes.end()) {
        NodeType* node = *it;
        if ((node->*filter)(coord))
            nodes.push_back(node);
        nof_visited ++;
        it ++;
    }
    // descend to children
    for (int i = 0; i < 8; i ++)
        if (octant->_children[i])
            find_nodes_by_point(coord, nodes, nof_visited, filter, octant->_children[i]);
}

template <typename NodeType>
NodeType* LooseOctree<NodeType>::find_one_node_by_point(const Vector3& coord,
                                                int& nof_visited,
                                                LooseOctree<NodeType>::filter filter,
                                                Octant<NodeType>* octant) {
    NodeType* node;
    
    // root is the default octant
    if (octant == 0)
        octant = _root;
    // return if coord is outside octant (including overlapping region)
    Vector3 dist = coord - octant->_center;
    double ext_length = 2.0 * octant->_length;
    if (fabs(dist.x) > ext_length ||
        fabs(dist.y) > ext_length ||
        fabs(dist.z) > ext_length)
        return 0;
    // check nodes stored in current Octant
    typename std::vector<NodeType *>::iterator it = octant->_nodes.begin();
    while (it != octant->_nodes.end()) {
        node = *it;
        nof_visited ++;
        it ++;
        if ((node->*filter)(coord))
            return node;
    }
    // descend to children
    for (int i = 0; i < 8; i ++) {
        if (octant->_children[i]) {
            node = find_one_node_by_point(coord, nof_visited, filter, octant->_children[i]);
            // if a node was found abort descend to other children
            if (node)
                return node;
        }
    }
    // Nothing found, return 0
    return 0;
}

template <typename NodeType>
void LooseOctree<NodeType>::get_statistics(int& max_depth,
                                 double& avg_node_depth,
                                 int& nof_octants,
                                 int& nof_nodes) {
    Octant<NodeType>* octant = _root;
    int depth = 0;
    double sum_depth = 0.0;

    max_depth = 0;
    nof_octants = 0;
    nof_nodes = 0;
    get_statistics_recursive(octant, depth, max_depth, sum_depth,
                             nof_octants, nof_nodes);
    avg_node_depth = sum_depth / nof_nodes;
}
  
template <typename NodeType>
void LooseOctree<NodeType>::get_statistics_recursive(Octant<NodeType>* octant,
                                           int depth,
                                           int& max_depth,
                                           double& sum_depth,
                                           int& nof_octants,
                                           int& nof_nodes) {
    // update statistics
    if (depth > max_depth)
        max_depth = depth;
    sum_depth += octant->_nodes.size() * depth;
    nof_octants ++;
    nof_nodes += octant->_nodes.size();

    // descend to children
    for (int i = 0; i < 8; i ++)
        if (octant->_children[i])
            get_statistics_recursive(octant->_children[i], depth+1, max_depth, sum_depth,
                                     nof_octants, nof_nodes);
}

template <typename NodeType>
void LooseOctree<NodeType>::export_vtk() {
    // get max depth
    int max_depth, nof_octants, nof_nodes;
    double avg_node_depth;
    get_statistics(max_depth, avg_node_depth, nof_octants, nof_nodes);

    // export all levels to separate VTK files
    for (int level = 0; level <= max_depth; ++ level) {
        std::ostringstream buf;
        buf << "octree_l" << level << ".vtk";
        export_vtk_level(level, buf.str());
    }
}

template <typename NodeType>
void LooseOctree<NodeType>::export_vtk_level(int level, std::string filename) {
    std::vector<int> cells;
    // recursive descend
    export_vtk_recursive(_root, 0, 0, 0, 0, level, cells);
    int nof_octants = cells.size() / 9;

    std::ofstream of;
    of.open(filename.c_str());
    rAssert(of.is_open());

    // VTK header
    of.precision(6);
    of << scientific;
    of << "# vtk DataFile Version 2.0" << endl;
    of << "generated using VtkExport::export_eigenfields" << endl;
    of << "ASCII" << endl << endl;
    of << "DATASET UNSTRUCTURED_GRID" << endl;

    // export all points in level, (2**depth + 1)**3 points
    int n = (1 << level) + 1;
    Vector3 origin = _root->_center;
    origin -= _root->_length;
    double side = 2.0*_root->_length / (1 << level);
    of << "POINTS " << n*n*n << " float" << endl;

    for (int z = 0; z < n; ++ z)
        for (int y = 0; y < n; ++ y)
            for (int x = 0; x < n; ++ x)
                of << origin.x + x*side << " "
                   << origin.y + y*side << " "
                   << origin.z + z*side << endl;
    of << endl;

    // cells
    of << "CELLS " << nof_octants << " " << cells.size() << endl;
    std::vector<int>::iterator it = cells.begin();
    for (int i = 0; it != cells.end(); ++ it, ++ i) {
        of << *it;
        if (i % 9 == 8)
            of << endl;
        else
            of << " ";
    }
    of << endl;

    // cell types
    of << "CELL_TYPES " << nof_octants << endl;
    for (int i = 0; i < nof_octants; ++ i)
        of << "11" << endl;
    of << endl;

    // cell data
    of << "CELL_DATA " << nof_octants << endl;
    of << "SCALARS cell_scalars int 1" << endl;
    of << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < nof_octants; ++ i)
        of << "0" << endl;
    of << endl;
}

template <typename NodeType>
void LooseOctree<NodeType>::export_vtk_recursive(Octant<NodeType>* octant,
                                       int depth,
                                       int x,
                                       int y,
                                       int z,
                                       int target_depth,
                                       std::vector<int>& cells) {
    if (depth == target_depth) {
        // export octant
        int n = (1 << depth) + 1;
        cells.push_back(8);
        cells.push_back((x+0) + (y+0)*n + (z+0)*n*n);
        cells.push_back((x+1) + (y+0)*n + (z+0)*n*n);
        cells.push_back((x+0) + (y+1)*n + (z+0)*n*n);
        cells.push_back((x+1) + (y+1)*n + (z+0)*n*n);
        cells.push_back((x+0) + (y+0)*n + (z+1)*n*n);
        cells.push_back((x+1) + (y+0)*n + (z+1)*n*n);
        cells.push_back((x+0) + (y+1)*n + (z+1)*n*n);
        cells.push_back((x+1) + (y+1)*n + (z+1)*n*n);
    } else {
        // descend to children
        for (int i = 0; i < 8; i ++) {
            if (octant->_children[i])
                export_vtk_recursive(octant->_children[i],
                                     depth+1,
                                     x + ((i & 0x1) == 0x1) * (1 << (target_depth - depth - 1)),
                                     y + ((i & 0x2) == 0x2) * (1 << (target_depth - depth - 1)),
                                     z + ((i & 0x4) == 0x4) * (1 << (target_depth - depth - 1)),
                                     target_depth,
                                     cells);
        }
    }
}

} // namespace mesh

#endif
