/***************************************************************************
                          octant.h  -  description
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

#ifndef OCTANT_H
#define OCTANT_H

#include <vector>
#include <ostream>
#include "vector3.h"

namespace mesh {

/** 
 * An octant in a LooseOctree
 * @author Roman Geus
 */

template <typename NodeType>
class Octant {
public:
    /** Constructor for the root Octant
     */
    Octant(const Vector3& center, double length);
    /** Constructor for non-root Octants
     */
    Octant(Octant* parent, int idx);
    /** Destructor: deallocates children Octants
     */
    ~Octant();
    /** Return child octant index corresponding to "node_center"
     */
    int get_child_index_for_point(const Vector3& node_center);
    /** Add node to this octant
     */
    void add_node(NodeType* node);
    /** Dump octant info including ancestor info
     */
    void dump(std::ostream& str);

    /** center of octant, x-, y- and z-coordinates */
    Vector3 _center;
    /** half of the side length of the octant */
    double _length;
    /** Parent node */
    Octant* _parent;
    /** 8 child octants */
    Octant* _children[8];
    /** Nodes stored in the octant */
    std::vector<NodeType *> _nodes;
};

template <typename NodeType>
Octant<NodeType>::Octant(const Vector3& center, double length)
    : _parent(0)
{
    // initialize all children to null.
    for (int i = 0; i < 8; i ++)
        _children[i] = 0;
    // set bounding box
    _center = center;
    _length = length;
}

template <typename NodeType>
Octant<NodeType>::Octant(Octant* parent, int idx)
    : _parent(parent)
{
    // initialize all children to null.
    for (int i = 0; i < 8; i ++)
        _children[i] = 0;
    // set bounding box
    _length = 0.5 * parent->_length;
    for (int i = 0; i < 3; i ++)
        if (idx & (0x1 << i))
            _center[i] = parent->_center[i] + _length;
        else
            _center[i] = parent->_center[i] - _length;
}

template <typename NodeType>
Octant<NodeType>::~Octant() {
    for (int i = 0; i < 8; i ++)
        delete _children[i];
    _parent = 0;  
}

template <typename NodeType>
int Octant<NodeType>::get_child_index_for_point(const Vector3& point) {
    int child_idx = 0;
    Vector3 disp = point - _center;
    if ( disp.x > 0 )
        child_idx += 1;
    if ( disp.y > 0 )
        child_idx += 2;
    if ( disp.z > 0 )
        child_idx += 4;
    return child_idx;
}

template <typename NodeType>
void Octant<NodeType>::add_node(NodeType* node) {
    _nodes.push_back(node);    
}

template <typename NodeType>
void Octant<NodeType>::dump(std::ostream& str) {
    Octant* octant = this;
        
    do {
        str << "Octant([" 
            << octant->_center.x << "," 
            << octant->_center.y << "," 
            << octant->_center.z << "," 
            << "], " << octant->_length << ")\n";
        octant = octant->_parent;
    } while (octant != 0);
}

} // namespace mesh

#endif
