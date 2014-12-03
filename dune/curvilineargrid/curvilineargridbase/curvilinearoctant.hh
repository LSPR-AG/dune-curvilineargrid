/***************************************************************************
                          curvilinearoctant.hh  -  description
                             -------------------
    begin                : Wed Feb 25 2004
    copyright            : (C) 2004 by Roman Geus
    email                : roman.geus@psi.ch

    edit                 : Mon Dec 1 2014
    author               : Aleksejs Fomins
    action               : Conversion to Dune Standard
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DUNE_CURVILINEAROCTANT_HH
#define DUNE_CURVILINEAROCTANT_HH

#include <vector>
#include <ostream>

#include <dune/common/fvector.hh>




namespace Dune {


template <class ct, int cdim, typename NodeType>
class CurvilinearOctant {
public:

	// Typedefs
	// ***************************************************

	typedef Dune::FieldVector<ct, cdim>   Vertex;



	// Construction
	// ***************************************************

    /** Constructor for the root CurvilinearOctant */
	CurvilinearOctant(const Vertex& center, double length) : parent_(0)
	{
	    // initialize all children to null.
	    for (int i = 0; i < 8; i++)  { children_[i] = 0; }

	    // set bounding box
	    center_ = center;
	    length_ = length;
	}


    /** Constructor for non-root Octants */
	CurvilinearOctant(CurvilinearOctant* parent, int childindex) : parent_(parent)
    {
        // initialize all children to null.
        for (int i = 0; i < 8; i ++)  { children_[i] = 0; }

        // set bounding box
        length_ = 0.5 * parent->length_;

        for (int i = 0; i < 3; i++)
        {
        	// childindex is 3 booleans stuck together into an integer.
        	// Each boolean for one of coordinates (x,y,z)
        	// Determines which of the possible 8 sub-octants this CurvilinearOctant is.
            if (childindex & (0x1 << i))  { center_[i] = parent->center_[i] + length_; }
            else                          { center_[i] = parent->center_[i] - length_; }
        }
    }


    /** Destructor: deallocates children Octants
     */
    ~CurvilinearOctant()
    {
        for (int i = 0; i < 8; i ++)  { delete children_[i]; }
        parent_ = 0;
    }



	// Implementation
	// ***************************************************

    /** Return child CurvilinearOctant index corresponding to "nodecenter_"
     */
    int childIndex(const Vertex& nodecenter_)
    {
        int rez = 0;
        Vertex disp = nodecenter_ - center_;

        if ( disp[0] > 0 )  { rez += 1; }  // x
        if ( disp[1] > 0 )  { rez += 2; }  // y
        if ( disp[2] > 0 )  { rez += 4; }  // z
        return rez;
    }


    /** Add node to this CurvilinearOctant
     */
    void addNode(NodeType* node)  { node_.push_back(node); }


    /** Dump CurvilinearOctant info including ancestor info
     */
    void dump(std::ostream& str)
    {
        CurvilinearOctant* CurvilinearOctant = this;

        do {
            str << "CurvilinearOctant(["
                << CurvilinearOctant->center_[0] << ","
                << CurvilinearOctant->center_[1] << ","
                << CurvilinearOctant->center_[2] << ","
                << "], " << CurvilinearOctant->length_ << ")\n";
            CurvilinearOctant = CurvilinearOctant->parent_;
        } while (CurvilinearOctant != 0);
    }






    /** center of CurvilinearOctant, x-, y- and z-coordinates */
    Vertex center_;
    /** half of the side length of the CurvilinearOctant */
    double length_;
    /** Parent node */
    CurvilinearOctant* parent_;
    /** 8 child octants */
    CurvilinearOctant* children_[8];
    /** Nodes stored in the CurvilinearOctant */
    std::vector<NodeType *> node_;
};


} // namespace Dune

#endif
