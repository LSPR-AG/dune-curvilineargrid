/** \file
 *  \brief Gmsh script for the standard sphere with radius = 1 length unit.
 *
 *  Copyright by Patrick Leidenberger and Benedikt Oswald, 2002-2009.
 *  All rights reserved.
 *
 *  Objective: Gmsh scriptum for the standard sphere with radius = 1 length unit. You can mesh with gmsh
 *  this file with: gmsh -3 sphere.geo
 *
 *  \author    Patrick Leidenberger, Benedikt Oswald
 *  \date      2006 dec 15, created, Patrick Leidenberger
 *  \date      2006 dec 15, modified, Patrick Leidenberger
 *  \date      2009 aug 03, modified, Benedikt Oswald, adapted parameters
 *  \date	   2009 sep 09, modified, Benedikt Oswald, changed to the standard sphere
 * 
 *  \warning   None.
 *  \attention 
 *  \bug
 *  \todo      
 */

scaling = 1.0;				// Scaling factor for all.
cl      = 0.2;				// Characteristic length for all.

// physical id's
abc1st = 101;
vacuum = 501;

// Center of sphere.
centerX = scaling * 0.0;
centerY = scaling * 0.0;
centerZ = scaling * 0.0;
 
// Radius of the sphere 
radius = 0.1;

// Create points on the sphere surface.
ipt1 = newp; Point(ipt1) = {centerX, centerY, centerZ, cl};
ipt2 = newp; Point(ipt2) = {centerX - radius, centerY, centerZ, cl};
ipt3 = newp; Point(ipt3) = {centerX + radius, centerY, centerZ, cl};
ipt4 = newp; Point(ipt4) = {centerX, centerY - radius, centerZ, cl};
ipt5 = newp; Point(ipt5) = {centerX, centerY + radius, centerZ, cl};
ipt6 = newp; Point(ipt6) = {centerX, centerY, centerZ - radius, cl};
ipt7 = newp; Point(ipt7) = {centerX, centerY, centerZ + radius, cl};

// Create circle sections connecting two points on sphere surface.
icl1  = newreg; Circle(icl1)  = {ipt2,ipt1,ipt4};
icl2  = newreg; Circle(icl2)  = {ipt2,ipt1,ipt5};
icl3  = newreg; Circle(icl3)  = {ipt2,ipt1,ipt6};
icl4  = newreg; Circle(icl4)  = {ipt2,ipt1,ipt7};
icl5  = newreg; Circle(icl5)  = {ipt3,ipt1,ipt4};
icl6  = newreg; Circle(icl6)  = {ipt3,ipt1,ipt5};
icl7  = newreg; Circle(icl7)  = {ipt3,ipt1,ipt6};
icl8  = newreg; Circle(icl8)  = {ipt3,ipt1,ipt7};
icl9  = newreg; Circle(icl9)  = {ipt4,ipt1,ipt6};
icl10 = newreg; Circle(icl10) = {ipt4,ipt1,ipt7};
icl11 = newreg; Circle(icl11) = {ipt5,ipt1,ipt6};
icl12 = newreg; Circle(icl12) = {ipt5,ipt1,ipt7};

// Make the surface mesh for 3 closed circle sections.
ill1 = newreg; Line Loop(ill1) = {icl1,-icl3,icl9} ; irs1= newreg; Ruled Surface(irs1) = {ill1};
ill2 = newreg; Line Loop(ill2) = {icl1,icl10,-icl4} ; irs2= newreg; Ruled Surface(irs2) = {ill2};
ill3 = newreg; Line Loop(ill3) = {icl2,-icl3,icl11} ; irs3= newreg; Ruled Surface(irs3) = {ill3};
ill4 = newreg; Line Loop(ill4) = {icl2,icl12,-icl4} ; irs4= newreg; Ruled Surface(irs4) = {ill4};
ill5 = newreg; Line Loop(ill5) = {icl5,icl9,-icl7}  ; irs5= newreg; Ruled Surface(irs5) = {ill5};
ill6 = newreg; Line Loop(ill6) = {icl5,icl10,-icl8} ; irs6= newreg; Ruled Surface(irs6) = {ill6};
ill7 = newreg; Line Loop(ill7) = {icl6,icl11,-icl7} ; irs7= newreg; Ruled Surface(irs7) = {ill7};
ill8 = newreg; Line Loop(ill8) = {icl6,icl12,-icl8} ; irs8= newreg; Ruled Surface(irs8) = {ill8};

ilsl1 = newreg;
Surface Loop(ilsl1) = {ill1+1,ill2+1,ill3+1,ill4+1,ill5+1,ill6+1,ill7+1,ill8+1};


// Define physcial surface
Physical Surface (abc1st) = {irs1,irs2,irs3,irs4,irs5,irs6,irs7,irs8};

// Define volume of sphere.
ivl1 = newv; Volume(ivl1) = {ilsl1};


// Define physical volume
Physical Volume ( vacuum ) = {ivl1};

