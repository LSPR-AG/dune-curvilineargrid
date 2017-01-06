
gridsize = 2.0; // prescribed mesh element size

rCube = {5.0, 7.0};
hCube = {0.0, 4.0, 7.0};

// Points LVL 1
p_lv1[0] = newp;  Point(p_lv1[0]) = {-rCube[0], -rCube[1], hCube[0], gridsize};
p_lv1[1] = newp;  Point(p_lv1[1]) = {rCube[0], -rCube[1], hCube[0], gridsize};
p_lv1[2] = newp;  Point(p_lv1[2]) = {rCube[0], rCube[1], hCube[0], gridsize};
p_lv1[3] = newp;  Point(p_lv1[3]) = {-rCube[0], rCube[1], hCube[0], gridsize};

// Points LVL 1
p_lv2[0] = newp;  Point(p_lv2[0]) = {-rCube[0], -rCube[1], hCube[1], gridsize / 3.0};
p_lv2[1] = newp;  Point(p_lv2[1]) = {rCube[0], -rCube[1], hCube[1], gridsize / 3.0};
p_lv2[2] = newp;  Point(p_lv2[2]) = {rCube[0], rCube[1], hCube[1], gridsize / 3.0};
p_lv2[3] = newp;  Point(p_lv2[3]) = {-rCube[0], rCube[1], hCube[1], gridsize / 3.0};

p_lv3[0] = newp;  Point(p_lv3[0]) = {-rCube[0], -rCube[1], hCube[2], gridsize};
p_lv3[1] = newp;  Point(p_lv3[1]) = {rCube[0], -rCube[1], hCube[2], gridsize};
p_lv3[2] = newp;  Point(p_lv3[2]) = {rCube[0], rCube[1], hCube[2], gridsize};
p_lv3[3] = newp;  Point(p_lv3[3]) = {-rCube[0], rCube[1], hCube[2], gridsize};

// Lines - Horizontal - LVL 1
linh_lv1[0] = newl;  Line(linh_lv1[0]) = {p_lv1[0], p_lv1[1]}; 
linh_lv1[1] = newl;  Line(linh_lv1[1]) = {p_lv1[1], p_lv1[2]};
linh_lv1[2] = newl;  Line(linh_lv1[2]) = {p_lv1[2], p_lv1[3]};
linh_lv1[3] = newl;  Line(linh_lv1[3]) = {p_lv1[3], p_lv1[0]};

// Lines - Horizontal - LVL 2
linh_lv2[0] = newl;  Line(linh_lv2[0]) = {p_lv2[0], p_lv2[1]}; 
linh_lv2[1] = newl;  Line(linh_lv2[1]) = {p_lv2[1], p_lv2[2]};
linh_lv2[2] = newl;  Line(linh_lv2[2]) = {p_lv2[2], p_lv2[3]};
linh_lv2[3] = newl;  Line(linh_lv2[3]) = {p_lv2[3], p_lv2[0]};

// Lines - Horizontal - LVL 3
linh_lv3[0] = newl;  Line(linh_lv3[0]) = {p_lv3[0], p_lv3[1]}; 
linh_lv3[1] = newl;  Line(linh_lv3[1]) = {p_lv3[1], p_lv3[2]};
linh_lv3[2] = newl;  Line(linh_lv3[2]) = {p_lv3[2], p_lv3[3]};
linh_lv3[3] = newl;  Line(linh_lv3[3]) = {p_lv3[3], p_lv3[0]};

// Lines - Vertical - LVL 1
linv_lv1[0] = newl;  Line(linv_lv1[0]) = {p_lv1[0], p_lv2[0]}; 
linv_lv1[1] = newl;  Line(linv_lv1[1]) = {p_lv1[1], p_lv2[1]};
linv_lv1[2] = newl;  Line(linv_lv1[2]) = {p_lv1[2], p_lv2[2]};
linv_lv1[3] = newl;  Line(linv_lv1[3]) = {p_lv1[3], p_lv2[3]};

// Lines - Vertical - LVL 2
linv_lv2[0] = newl;  Line(linv_lv2[0]) = {p_lv2[0], p_lv3[0]}; 
linv_lv2[1] = newl;  Line(linv_lv2[1]) = {p_lv2[1], p_lv3[1]};
linv_lv2[2] = newl;  Line(linv_lv2[2]) = {p_lv2[2], p_lv3[2]};
linv_lv2[3] = newl;  Line(linv_lv2[3]) = {p_lv2[3], p_lv3[3]};

Periodic Line {linh_lv1[0]} = {linh_lv3[0]};
Periodic Line {linh_lv1[1]} = {linh_lv3[1]};
Periodic Line {linh_lv1[2]} = {linh_lv3[2]};
Periodic Line {linh_lv1[3]} = {linh_lv3[3]};

Periodic Line {linh_lv1[0]} = {-linh_lv1[2]};
Periodic Line {linh_lv1[1]} = {-linh_lv1[3]};
Periodic Line {linh_lv2[0]} = {-linh_lv2[2]};
Periodic Line {linh_lv2[1]} = {-linh_lv2[3]};
Periodic Line {linh_lv3[0]} = {-linh_lv3[2]};
Periodic Line {linh_lv3[1]} = {-linh_lv3[3]};

Periodic Line {linv_lv1[0]} = {linv_lv1[3]};
Periodic Line {linv_lv1[1]} = {linv_lv1[2]};
Periodic Line {linv_lv1[1]} = {linv_lv1[0]};
Periodic Line {linv_lv1[2]} = {linv_lv1[3]};

Periodic Line {linv_lv2[0]} = {linv_lv2[3]};
Periodic Line {linv_lv2[1]} = {linv_lv2[2]};
Periodic Line {linv_lv2[1]} = {linv_lv2[0]};
Periodic Line {linv_lv2[2]} = {linv_lv2[3]};

llh_lv1[0] = newreg;  Line Loop(llh_lv1[0]) = {linh_lv1[]};
llh_lv2[0] = newreg;  Line Loop(llh_lv2[0]) = {linh_lv2[]};
llh_lv3[0] = newreg;  Line Loop(llh_lv3[0]) = {linh_lv3[]};

llv_lv1[0] = newreg;  Line Loop(llv_lv1[0]) = {linh_lv1[0], linv_lv1[1], -linh_lv2[0], -linv_lv1[0]};
llv_lv1[1] = newreg;  Line Loop(llv_lv1[1]) = {linh_lv1[1], linv_lv1[2], -linh_lv2[1], -linv_lv1[1]};
llv_lv1[2] = newreg;  Line Loop(llv_lv1[2]) = {linh_lv1[2], linv_lv1[3], -linh_lv2[2], -linv_lv1[2]};
llv_lv1[3] = newreg;  Line Loop(llv_lv1[3]) = {linh_lv1[3], linv_lv1[0], -linh_lv2[3], -linv_lv1[3]};

llv_lv2[0] = newreg;  Line Loop(llv_lv2[0]) = {linh_lv2[0], linv_lv2[1], -linh_lv3[0], -linv_lv2[0]};
llv_lv2[1] = newreg;  Line Loop(llv_lv2[1]) = {linh_lv2[1], linv_lv2[2], -linh_lv3[1], -linv_lv2[1]};
llv_lv2[2] = newreg;  Line Loop(llv_lv2[2]) = {linh_lv2[2], linv_lv2[3], -linh_lv3[2], -linv_lv2[2]};
llv_lv2[3] = newreg;  Line Loop(llv_lv2[3]) = {linh_lv2[3], linv_lv2[0], -linh_lv3[3], -linv_lv2[3]};

sh_lv1[0] = news;  Plane Surface(sh_lv1[0]) = {llh_lv1[0]};
sh_lv2[0] = news;  Plane Surface(sh_lv2[0]) = {llh_lv2[0]};
sh_lv3[0] = news;  Plane Surface(sh_lv3[0]) = {llh_lv3[0]};
sv_lv1[0] = news;  Plane Surface(sv_lv1[0]) = {llv_lv1[0]};
sv_lv1[1] = news;  Plane Surface(sv_lv1[1]) = {llv_lv1[1]};
sv_lv1[2] = news;  Plane Surface(sv_lv1[2]) = {llv_lv1[2]};
sv_lv1[3] = news;  Plane Surface(sv_lv1[3]) = {llv_lv1[3]};
sv_lv2[0] = news;  Plane Surface(sv_lv2[0]) = {llv_lv2[0]};
sv_lv2[1] = news;  Plane Surface(sv_lv2[1]) = {llv_lv2[1]};
sv_lv2[2] = news;  Plane Surface(sv_lv2[2]) = {llv_lv2[2]};
sv_lv2[3] = news;  Plane Surface(sv_lv2[3]) = {llv_lv2[3]};

Periodic Surface(sh_lv1[0]) {linh_lv1[]} = sh_lv3[0] {linh_lv3[]};
Periodic Surface(sv_lv1[0]) {linh_lv1[0], linv_lv1[1], -linh_lv2[0], -linv_lv1[0]} = sv_lv1[2] {-linh_lv1[2], linv_lv1[2], linh_lv2[2], -linv_lv1[3] };
Periodic Surface(sv_lv1[1]) {linh_lv1[1], linv_lv1[2], -linh_lv2[1], -linv_lv1[1]} = sv_lv1[3] {-linh_lv1[3], linv_lv1[3], linh_lv2[3], -linv_lv1[0] };
Periodic Surface(sv_lv2[0]) {linh_lv2[0], linv_lv2[1], -linh_lv3[0], -linv_lv2[0]} = sv_lv2[2] {-linh_lv2[2], linv_lv2[2], linh_lv3[2], -linv_lv2[3] };
Periodic Surface(sv_lv2[1]) {linh_lv2[1], linv_lv2[2], -linh_lv3[1], -linv_lv2[1]} = sv_lv2[3] {-linh_lv2[3], linv_lv2[3], linh_lv3[3], -linv_lv2[0] };

//Periodic Surface {sv_lv1[0]} = {-sv_lv1[2]};
//Periodic Surface {sv_lv1[1]} = {-sv_lv1[3]};

sl[0] = newreg;  Surface Loop(sl[0]) = {sh_lv1[], sh_lv2[], sv_lv1[]};
sl[1] = newreg;  Surface Loop(sl[1]) = {sh_lv2[], sh_lv3[], sv_lv2[]};
vol[0] = newv;   Volume(vol[0]) = {sl[0]};
vol[1] = newv;   Volume(vol[1]) = {sl[1]};


// Define Tags
abc1st = 10000;
box1 = 1000;
box2 = 2000;

// Define physcial surface
Physical Surface (abc1st) = {sh_lv1[], sh_lv3[], sv_lv1[], sv_lv2[]};

// Define physical volume
Physical Volume (box1) = {vol[0]};
Physical Volume (box2) = {vol[1]};



