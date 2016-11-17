/** \file
 *  \brief Gmsh script that generates three two spheres with outermost radius=1 unit
 *
 *  Copyright by Patrick Leidenberger and Benedikt Oswald, 2006-2009.
 *  All rights reserved.
 *
 *  Objective: Gmsh script of a sphere in a sphere, the classical mesh for calculating
 *             the scattering of an incoming plane wave at a dielectric or metallic sphere,
 *             aka. Mie scattering.
 *
 *  \author    Patrick Leidenberger, Benedikt Oswald
 *  \date      2006 dec 15, created, Patrick Leidenberger.
 *  \date      2006 dec 15, modified, Patrick Leidenberger,
 *  \date      2007 apr 26, modified, Benedikt Oswald, added physical boundary and volume tags
 *  \date      2009 aug 04, modified, Benedikt Oswald, changes to boundary and volume tags
 *  \date      2009 aug 04, renamed, Benedikt Oswald, adapted to two nested spheres.
 * 
 *  \warning   None.
 *  \attention 
 *  \bug
 *  \todo      
 */

// define physical surface and volume tags
surf_in     = 101;
surf_out    = 102;
surf_db     = 103;
vol_in      = 501;
vol_shell   = 502;
vol_out     = 503;


// Radius of the sphere 
radius0 = 0.05;				// radius of inner most sphere
radius1 = 0.50;				// outer radius of surrounding shell
radius2 = 1.00;                         // radius of outer most sphere

cl0     = radius0 * 0.50;		// Characteristic length sphere 0 = inner most sphere
cl1     = radius1 * 0.50;		// Characteristic length sphere 1 = sphere in between inner most and outer most sphere
cl2     = radius2 * 0.50;		// Characteristic length sphere 2 = outer most sphere

// Center of sphere.
centerX = 0.0;
centerY = 0.0;
centerZ = 0.0;

// Create points on sphere 0 surface = innermost sphere
ipt1  = newp; Point(ipt1)  = {centerX, centerY, centerZ, cl0};
ipt2  = newp; Point(ipt2)  = {centerX - radius0, centerY, centerZ, cl0};
ipt3  = newp; Point(ipt3)  = {centerX + radius0, centerY, centerZ, cl0};
ipt4  = newp; Point(ipt4)  = {centerX, centerY - radius0, centerZ, cl0};
ipt5  = newp; Point(ipt5)  = {centerX, centerY + radius0, centerZ, cl0};
ipt6  = newp; Point(ipt6)  = {centerX, centerY, centerZ - radius0, cl0};
ipt7  = newp; Point(ipt7)  = {centerX, centerY, centerZ + radius0, cl0};

// Create points on sphere 1 surface = sphere in between
ipt8  = newp; Point(ipt8)  = {centerX, centerY, centerZ, cl1};
ipt9  = newp; Point(ipt9)  = {centerX - radius1, centerY, centerZ, cl1};
ipt10 = newp; Point(ipt10) = {centerX + radius1, centerY, centerZ, cl1};
ipt11 = newp; Point(ipt11) = {centerX, centerY - radius1, centerZ, cl1};
ipt12 = newp; Point(ipt12) = {centerX, centerY + radius1, centerZ, cl1};
ipt13 = newp; Point(ipt13) = {centerX, centerY, centerZ - radius1, cl1};
ipt14 = newp; Point(ipt14) = {centerX, centerY, centerZ + radius1, cl1};

// Create points on sphere 2 surface = outermost sphere
ipt15 = newp; Point(ipt15) = {centerX, centerY, centerZ, cl2};
ipt16 = newp; Point(ipt16) = {centerX - radius2, centerY, centerZ, cl2};
ipt17 = newp; Point(ipt17) = {centerX + radius2, centerY, centerZ, cl2};
ipt18 = newp; Point(ipt18) = {centerX, centerY - radius2, centerZ, cl2};
ipt19 = newp; Point(ipt19) = {centerX, centerY + radius2, centerZ, cl2};
ipt20 = newp; Point(ipt20) = {centerX, centerY, centerZ - radius2, cl2};
ipt21 = newp; Point(ipt21) = {centerX, centerY, centerZ + radius2, cl2};

// Create circle sections connecting two points on sphere 0 surface.
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
// Create circle sections connecting two points on sphere 1 surface.
icl13 = newreg; Circle(icl13) = {ipt9,ipt8,ipt11};
icl14 = newreg; Circle(icl14) = {ipt9,ipt8,ipt12};
icl15 = newreg; Circle(icl15) = {ipt9,ipt8,ipt13};
icl16 = newreg; Circle(icl16) = {ipt9,ipt8,ipt14};
icl17 = newreg; Circle(icl17) = {ipt10,ipt8,ipt11};
icl18 = newreg; Circle(icl18) = {ipt10,ipt8,ipt12};
icl19 = newreg; Circle(icl19) = {ipt10,ipt8,ipt13};
icl20 = newreg; Circle(icl20) = {ipt10,ipt8,ipt14};
icl21 = newreg; Circle(icl21) = {ipt11,ipt8,ipt13};
icl22 = newreg; Circle(icl22) = {ipt11,ipt8,ipt14};
icl23 = newreg; Circle(icl23) = {ipt12,ipt8,ipt13};
icl24 = newreg; Circle(icl24) = {ipt12,ipt8,ipt14};
// Create circle sections connecting two points on sphere 2 surface.
icl25 = newreg; Circle(icl25) = {ipt16,ipt15,ipt18};
icl26 = newreg; Circle(icl26) = {ipt16,ipt15,ipt19};
icl27 = newreg; Circle(icl27) = {ipt16,ipt15,ipt20};
icl28 = newreg; Circle(icl28) = {ipt16,ipt15,ipt21};
icl29 = newreg; Circle(icl29) = {ipt17,ipt15,ipt18};
icl30 = newreg; Circle(icl30) = {ipt17,ipt15,ipt19};
icl31 = newreg; Circle(icl31) = {ipt17,ipt15,ipt20};
icl32 = newreg; Circle(icl32) = {ipt17,ipt15,ipt21};
icl33 = newreg; Circle(icl33) = {ipt18,ipt15,ipt20};
icl34 = newreg; Circle(icl34) = {ipt18,ipt15,ipt21};
icl35 = newreg; Circle(icl35) = {ipt19,ipt15,ipt20};
icl36 = newreg; Circle(icl36) = {ipt19,ipt15,ipt21};

// Make the surface mesh for 3 closed circle sections; sphere0.
ill1 = newreg; Line Loop(ill1) = {icl1,-icl3,icl9} ; irs1= newreg; Ruled Surface(irs1) = {ill1};
ill2 = newreg; Line Loop(ill2) = {icl1,icl10,-icl4} ; irs2= newreg; Ruled Surface(irs2) = {ill2};
ill3 = newreg; Line Loop(ill3) = {icl2,-icl3,icl11} ; irs3= newreg; Ruled Surface(irs3) = {ill3};
ill4 = newreg; Line Loop(ill4) = {icl2,icl12,-icl4} ; irs4= newreg; Ruled Surface(irs4) = {ill4};
ill5 = newreg; Line Loop(ill5) = {icl5,icl9,-icl7} ; irs5= newreg; Ruled Surface(irs5) = {ill5};
ill6 = newreg; Line Loop(ill6) = {icl5,icl10,-icl8} ; irs6= newreg; Ruled Surface(irs6) = {ill6};
ill7 = newreg; Line Loop(ill7) = {icl6,icl11,-icl7} ; irs7= newreg; Ruled Surface(irs7) = {ill7};
ill8 = newreg; Line Loop(ill8) = {icl6,icl12,-icl8} ; irs8= newreg; Ruled Surface(irs8) = {ill8};
// Make the surface mesh for 3 closed circle sections; sphere1.
ill9 = newreg; Line Loop(ill9) = {icl13,-icl15,icl21} ; irs9= newreg; Ruled Surface(irs9) = {ill9};
ill10 = newreg; Line Loop(ill10) = {icl13,icl22,-icl16} ; irs10= newreg; Ruled Surface(irs10) = {ill10};
ill11 = newreg; Line Loop(ill11) = {icl14,-icl15,icl23} ; irs11= newreg; Ruled Surface(irs11) = {ill11};
ill12 = newreg; Line Loop(ill12) = {icl14,icl24,-icl16} ; irs12= newreg; Ruled Surface(irs12) = {ill12};
ill13 = newreg; Line Loop(ill13) = {icl17,icl21,-icl19} ; irs13= newreg; Ruled Surface(irs13) = {ill13};
ill14 = newreg; Line Loop(ill14) = {icl17,icl22,-icl20} ; irs14= newreg; Ruled Surface(irs14) = {ill14};
ill15 = newreg; Line Loop(ill15) = {icl18,icl23,-icl19} ; irs15= newreg; Ruled Surface(irs15) = {ill15};
ill16 = newreg; Line Loop(ill16) = {icl18,icl24,-icl20} ; irs16= newreg; Ruled Surface(irs16) = {ill16};
// Make the surface mesh for 3 closed circle sections; sphere2.
ill17 = newreg; Line Loop(ill17) = {icl25,-icl27,icl33} ; irs17= newreg; Ruled Surface(irs17) = {ill17};
ill18 = newreg; Line Loop(ill18) = {icl25,icl34,-icl28} ; irs18= newreg; Ruled Surface(irs18) = {ill18};
ill19 = newreg; Line Loop(ill19) = {icl26,-icl27,icl35} ; irs19= newreg; Ruled Surface(irs19) = {ill19};
ill20 = newreg; Line Loop(ill20) = {icl26,icl36,-icl28} ; irs20= newreg; Ruled Surface(irs20) = {ill20};
ill21 = newreg; Line Loop(ill21) = {icl29,icl33,-icl31} ; irs21= newreg; Ruled Surface(irs21) = {ill21};
ill22 = newreg; Line Loop(ill22) = {icl29,icl34,-icl32} ; irs22= newreg; Ruled Surface(irs22) = {ill22};
ill23 = newreg; Line Loop(ill23) = {icl30,icl35,-icl31} ; irs23= newreg; Ruled Surface(irs23) = {ill23};
ill24 = newreg; Line Loop(ill24) = {icl30,icl36,-icl32} ; irs24= newreg; Ruled Surface(irs24) = {ill24};

ilsl1 = newreg;
Surface Loop(ilsl1) = {ill1+1,ill2+1,ill3+1,ill4+1,ill5+1,ill6+1,ill7+1,ill8+1};
ilsl2 = newreg;
Surface Loop(ilsl2) = {ill1+1,ill2+1,ill3+1,ill4+1,ill5+1,ill6+1,ill7+1,ill8+1,ill9+1,ill10+1,ill11+1,ill12+1,ill13+1,ill14+1,ill15+1,ill16+1};
ilsl3 = newreg;
Surface Loop(ilsl3) = {ill9+1,ill10+1,ill11+1,ill12+1,ill13+1,ill14+1,ill15+1,ill16+1,ill17+1,ill18+1,ill19+1,ill20+1,ill21+1,ill22+1,ill23+1,ill24+1};


// define physical surface tags

Physical Surface (surf_in) = {irs1,irs2,irs3,irs4,irs5,irs6,irs7,irs8};
Physical Surface (surf_out) = {irs17,irs18,irs19,irs20,irs21,irs22,irs23,irs24};
Physical Surface (surf_db) = {irs9,irs10,irs11,irs12,irs13,irs14,irs15,irs16};


// define volume of sphere.
ivl1 = newv; Volume(ivl1) = {ilsl1};
ivl2 = newv; Volume(ivl2) = {ilsl2};
ivl3 = newv; Volume(ivl3) = {ilsl3};


// define physical volumes
Physical Volume (vol_in) = ivl1;
Physical Volume (vol_shell) = ivl2;
Physical Volume (vol_out) = ivl3;


