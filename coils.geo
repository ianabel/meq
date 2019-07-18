/*********************************************************************
 *
 *  Generate a rectangular space, then mesh it
 *
 *********************************************************************/

h = 2e-1;
coil_h = 1e-2;

Point(1) = {0.05, -1.5, 0, h};
Point(2) = {1.5, -1.5,  0, h} ;
Point(3) = {1.5, 1.5, 0, h} ;
Point(4) = {0.05, 1.5, 0, h} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Curve Loop(1) = {1,2,3,4};

Point(5) = {0.975,0.975,0,coil_h};
Point(6) = {1.025,0.975,0,coil_h};
Point(7) = {1.025,1.025,0,coil_h};
Point(8) = {0.975,1.025,0,coil_h};

Line(5) = {5,6} ;
Line(6) = {6,7} ;
Line(7) = {7,8} ;
Line(8) = {8,5} ;

Curve Loop(2) = {5,6,7,8};

Point(9) = {0.975,-0.975,0,coil_h};
Point(10) = {1.025,-0.975,0,coil_h};
Point(11) = {1.025,-1.025,0,coil_h};
Point(12) = {0.975,-1.025,0,coil_h};

Line(9) = {9,10} ;
Line(10) = {10,11} ;
Line(11) = {11,12} ;
Line(12) = {12,9} ;

Curve Loop(3) = {9,10,11,12};

Plane Surface(1) = {1,2,3};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Physical Surface("Domain",1) = {1,2,3};

