// Copyright (C) 2025 Andreas Stahel
// 
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see
// <https://www.gnu.org/licenses/>.


// Author: Andreas Stahel <andreas.stahel@gmx.com>
// Created: 2025-01-06

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/ov-struct.h>

#define NEUMANN   -2
#define DIRICHLET -1

DEFUN_DLD (FEMEquationCubic, args, ,
   "-*- texinfo -*-\n\
@deftypefn {} {} [@var{A},@var{b}] = FEMEquationCubic(@var{mesh},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2})\n\
\n\
sets up the system of linear equations for a numerical solution of a PDE\n\
using a triangular mesh with elements of order 3\n\
\n\
@verbatim\n\
-div(a*grad u - u*(bx,by))+ b0*u = f         in domain\n\
                               u = gD        on Dirichlet boundary\n\
        n*(a*grad u - u*(bx,by)) = gN1+g2N*u on Neumann boundary\n\
@end verbatim\n\
\n\
parameters:\n\
@itemize\n\
@item @var{mesh} triangular mesh of order 3 describing the domain and the boundary types\n\
@item @var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2}\n\
are the coefficients and functions describing the PDE.\n\
@*Any constant function can be given by its scalar value.\n\
@*The functions @var{a},@var{b0},@var{bx},@var{by} and @var{f} may also be given as vectors\n\
with the values of the function at the Gauss points.\n\
@item The coefficient @var{a} can als be a symmetric matrix @var{a=[axx,axy;axy,ayy]} given by the row vector @var{[axx,ayy,axy]}.\n\
@* It can be given as row vector or as string with the function name or as nx3 matrix with the values at the Gauss points.\
@end itemize\n\
\n\
return values:\n\
@itemize\n\
@item@var{A}, @var{b}: matrix and vector for the linear system to be solved, @var{A}*u-@var{b}=0\n\
@end itemize\n\
@c BEGIN_CUT_TEXINFO\n\
@seealso{FEMEquation, FEMEquationQuad, FEMEquationComplex, FEMEquationQuadComplex, FEMEquationCubicComplex, BVP2D, BVP2Dsym, BVP2eig, IBVP2D, I2BVP2D, CreateMeshRect, CreateMeshTriangle}\n\
@c END_CUT_TEXINFO\n\
@end deftypefn")

{
  octave_value_list retval;
  int nargin = args.length ();
  if (nargin !=9 ) {
    print_usage ();
    return retval;
  }
  using namespace std;
  octave_value_list argin;
  octave_value_list res;
  int nargout;
  
  /* INPUT LOADING */
  octave_map arg0 = args(0).map_value ();
  octave_map::const_iterator p1 = arg0.seek ("nodes");
  const Matrix nodes  =  arg0.contents(p1)(0).matrix_value();
  p1 = arg0.seek ("nodesT");
  const ColumnVector nodesT    =  arg0.contents(p1)(0).column_vector_value();
  p1 = arg0.seek ("edges");
  const Matrix edges           =  arg0.contents(p1)(0).matrix_value();
  p1 = arg0.seek ("edgesT");
  const ColumnVector edgesT    =  arg0.contents(p1)(0).column_vector_value();
  p1 = arg0.seek ("elem");
  const Matrix elem            =  arg0.contents(p1)(0).matrix_value();
  p1 = arg0.seek ("elemT");
  const ColumnVector elemT     =  arg0.contents(p1)(0).column_vector_value();
  p1 = arg0.seek ("elemArea");
  const ColumnVector elemArea  =  arg0.contents(p1)(0).column_vector_value();
  p1 = arg0.seek ("node2DOF");
  const ColumnVector n2d       =  arg0.contents(p1)(0).column_vector_value();
  p1 = arg0.seek ("nDOF");
  const int nDOF               =  arg0.contents(p1)(0).int_value();
  p1 = arg0.seek ("GP");
  const Matrix GP              =  arg0.contents(p1)(0).matrix_value();
  p1 = arg0.seek ("GPT");
  const ColumnVector GPT       =  arg0.contents(p1)(0).column_vector_value();
  int nGP = GP.rows();

  /* INPUT CHECK and Evaluate*/
  string Func;
  bool isotropic = true;

  //  evaluate function a
  ColumnVector aV, aVtmp;
  int ii;
  if (args(1).is_string()) {  // function given as string
    Func = args(1).string_value();
    argin(0) = GP;
    nargout  = 1;
    res = octave::feval (Func, argin, nargout);
    aV = res(0).column_vector_value();
    if (res(0).columns()==3)
      { isotropic = false;} // store all coefficients in one vector aV
                            // first a11, then a22, then a12
  }
  else if(args(1).is_real_scalar()){ //function given as scalar
    aV.resize(nGP);
    aV.fill(args(1).double_value());
  }
  else {// function given by its values
    if(args(1).is_real_scalar()){ //function given as scalar
      aV.resize(nGP);
      aV.fill(args(1).double_value());
    }
    else {
      //      octave_stdout <<"not a scalar\n";
      if (args(1).columns()==3)  // anisotropic
      { isotropic = false; // store all coefficients in one vector aV
                           // first a11, then a22, then a12
	if (args(1).rows()==1){// only three values given
	  aVtmp = args(1).column_vector_value();
	  aV.resize(3*nGP);
	  for(ii = 0; ii < nGP; ii = ii+1){
	    aV(ii)             = aVtmp(0);
	    aV(ii+nGP)   = aVtmp(1);
	    aV(ii+2*nGP) = aVtmp(2);
	  }
	}
	else {//  all a values in three columns
	  aV = args(1).column_vector_value();
	}
      }
      else {//  all a values in three columns
	aV = args(1).column_vector_value();
      }
    }
  }

  // evaluate function b
  ColumnVector bV; 
  if (args(2).is_string()) {  // function given as string
    Func = args(2).string_value();
    nargout  = 1;
    argin(0) = GP;
    res = octave::feval (Func, argin, nargout);
    bV = res(0).column_vector_value();
  }
  else if(args(2).is_real_scalar()){ //function given as scalar
    bV.resize(nGP);
    bV.fill(args(2).double_value());
  }
  else {  // function given by its values
    bV = args(2).column_vector_value();
  }
  
  bool convectionFlag = false;
  // evaluate function bx
  ColumnVector bxV; 
  if (args(3).is_string()) {  // function given as string
    convectionFlag = true;
    Func    = args(3).string_value();
    nargout = 1;
    argin(0) = GP;
    res = octave::feval (Func, argin, nargout);
    bxV=res(0).column_vector_value();
  }
  else if(args(3).is_real_scalar()){ //function given as scalar
    bxV.resize(nGP);
    bxV.fill(args(3).double_value());
    if (args(3).double_value()!=0){convectionFlag=true;}
  }
  else {  // function given by its values
    bxV = args(3).column_vector_value();
    convectionFlag=true;
 }

  // evaluate function by
  ColumnVector byV; 
  if (args(4).is_string()) {  // function given as string
    Func = args(4).string_value();
    convectionFlag = true;
    nargout  = 1;
    argin(0) = GP;
    res = octave::feval (Func, argin, nargout);
    byV = res(0).column_vector_value();
  }
  else if(args(4).is_real_scalar()){ //function given as scalar
    byV.resize(nGP);
    byV.fill(args(4).double_value());
    if (args(4).double_value() !=0){convectionFlag=true;}
  }
  else {  // function given by its values
    byV = args(4).column_vector_value();
    convectionFlag = true;
  }

  //  evaluate function f
  ColumnVector fV; 
  if (args(5).is_string()) {  // function given as string
    Func = args(5).string_value();
    nargout  = 1;
    argin(0) = GP;
    res = octave::feval (Func, argin, nargout);
    fV = res(0).column_vector_value();
  }
  else if(args(5).is_real_scalar()){ //function given as scalar
    fV.resize(nGP);
    fV.fill(args(5).double_value());
  }
  else {  // function given by its values
    fV = args(5).column_vector_value();
  }
  
//  string gDFunc  = args(6).string_value();
  string gDFunc;
  double gDValue = 0.0;;
  bool gDscalar  = true;
  if (args(6).is_string()) {  // function given as string
    gDFunc   = args(6).string_value();
    gDscalar = false;
  }
  else {gDValue = args(6).double_value();}

//  string gN1Func  = args(7).string_value();
  string gN1Func;
  double gN1Value = 0.0;
  bool gN1scalar  = true;
  if (args(7).is_string()) {  // function given as string
    gN1Func   = args(7).string_value();
    gN1scalar = false;
  }
  else {gN1Value = args(7).double_value();}

//  string gN2Func  = args(8).string_value();
  string gN2Func;
  double gN2Value = 0.0;
  bool gN2scalar  = true;
  if (args(8).is_string()) {  // function given as string
    gN2Func   = args(8).string_value();
    gN2scalar = false;
  }
  else {gN2Value = args(8).double_value();}
  
  /* END INIT */
  
  /* MEMORY ALLOC AND CONST */ 
  int ptrDOF = 0;
  
  // interpolation matrix for second and third order elements
  double l1 = (12.0-2.0*sqrt(15.0))/21.0;
  double l2 = (12.0+2.0*sqrt(15.0))/21.0;
  double w1 = (155.0-sqrt(15.0))/2400.0;
  double w2 = (155.0+sqrt(15.0))/2400.0;
  double w3 = 0.1125;

  ColumnVector w(7);  // weights of Gauss points
  w(0) = w1; w(1) = w1; w(2) = w1; w(3) = w2; w(4) = w2; w(5) = w2; w(6) = w3;
  ColumnVector xi(7); // first coordinates of Gauss points
  xi(0)=0.5*l1; xi(1)=1.0-l1; xi(2)=0.5*l1; xi(3)=0.5*l2; xi(4)=1.0-l2; xi(5)=0.5*l2; xi(6)=1.0/3.0;
  ColumnVector nu(7); // second coordinates of Gauss points
  nu(0)=0.5*l1; nu(1)=0.5*l1; nu(2)=1.0-l1; nu(3)=0.5*l2; nu(4)=0.5*l2; nu(5)=1.0-l2; nu(6)=1.0/3.0;

  double xit, nut, xinu;
  
  Matrix M(7,10);   // interpolate values of function 
  Matrix Mxi(7,10); // interpolate values of partial derivative w.r. to xi
  Matrix Mnu(7,10); // interpolate values of partial derivative w.r. to nu
  for (int row=0;row<7;row++){
    xit = xi(row); nut = nu(row); xinu = xit+nut;
    M(row,0) = (1.0-xinu)*(1.0-3.0*xinu)*(1.0-1.5*xinu); 
    M(row,1) = xit*(3.0*xit-1.0)*(1.5*xit-1.0);
    M(row,2) = nut*(3.0*nut-1.0)*(1.5*nut-1.0);
    M(row,3) = 4.5*xit*nut*(3.0*xit-1.0);
    M(row,4) = 4.5*xit*nut*(3.0*nut-1.0);
    M(row,5) = 4.5*nut*(1.0-xinu)*(3.0*nut-1.0);
    M(row,6) = 9.0*nut*(1.0-xinu)*(1.0-1.5*xinu);
    M(row,7) = 9.0*xit*(1.0-1.5*xinu)*(1.0-xinu);
    M(row,8) = 4.5*xit*(3.0*xit-1.0)*(1.0-xinu);
    M(row,9) = 27.0*xit*nut*(1.0-xinu);

    Mxi(row,0) = -5.5+18.0*xinu-13.5*xinu*xinu; 
    Mxi(row,1) = 1.0 - xit*(9.0-13.5*xit);
    Mxi(row,2) = 0.0;
    Mxi(row,3) = nut*(-4.5+27.0*xit);
    Mxi(row,4) = nut*(-4.5+13.5*nut);
    Mxi(row,5) = -Mxi(row,4);
    Mxi(row,6) = -22.5*nut + 27.0*nut*xinu;
    Mxi(row,7) = 9.0-45.0*xit- 22.5*nut +40.5*xit*xit+(54.0*xit+13.5*nut)*nut;
    Mxi(row,8) = -4.5+36.0*xit+4.5*nut - xit*(40.5*xit +27.0*nut);
    Mxi(row,9) = 27.0*(1.0- 2.0*xit - nut)*nut;
    
    Mnu(row,0) = Mxi(row,0);
    Mnu(row,1) = 0.0;
    Mnu(row,2) = 1.0 - 9.0*(1.0-1.5*nut)*nut;
    Mnu(row,3) = -4.5*(1.0-3.0*xit)*xit;
    Mnu(row,4) = -4.5*(1.0-6.0*nut)*xit;
    Mnu(row,5) = -4.5*(1.0-xit) + 36.0*nut -27.0*(xit+1.5*nut)*nut;
    Mnu(row,6) = 9.0-22.5*xit-45.0*nut + 13.5*xit*xit +(54.0*xit+40.5*nut)*nut;
    Mnu(row,7) = -22.5*xit + 27.0*xit*xinu;
    Mnu(row,8) = 4.5*xit*(1.0-3.0*xit);
    Mnu(row,9) = 27.0*xit*(1.0 - xit - 2.0*nut);
  };
  Matrix Mtrans = M.transpose();

  // diagonal matrices (some will be overwritten later, but have the same size as diagb, therefore they are initialized as follows)
  DiagMatrix diagb(7,7);
  DiagMatrix diagaxx(7,7);
  DiagMatrix diagayy(7,7);
  DiagMatrix diagaxy(7,7);
  DiagMatrix diagb1(7,7);
  DiagMatrix diagb2(7,7);
  // weights for edge contribution	
  ColumnVector wc(3);
  wc(0) = 5.0/18.0; wc(1) = 8.0/18.0;  wc(2) = 5.0/18.0;
  DiagMatrix diagwc(wc);

  // matrices and variables for element stiffness matrix
  Matrix Ab1(10,10);
  Matrix Ab2(10,10);
  Matrix corners(10,2);
  Matrix T(2,2);
  double detT;
  Matrix elMat (10,10);
  ColumnVector elVec (10);
  ColumnVector wvalVec (7);
  Matrix ttx(7,10);
  Matrix tty(7,10);
  Matrix Axx(10,10);
  Matrix Ayy(10,10);
  Matrix Axy(10,10);

  // for assembling of global stiffnes matrix
  RowVector tmpCorner(2); // for jus one corner
  ColumnVector dofs(10);  // degrees of freedom in element
  ColumnVector dofsB(4); // degrees of freedom on edge
  
  // matrices and variables for edge contribution
  Matrix MB1(4,4);
    MB1(0,0) = -1.0; MB1(0,1) = 9.0; MB1(0,2) = 9.0; MB1(0,3) = -1.0;
    MB1(1,0) = 2.0; MB1(1,1) = -54.0; MB1(1,2) = 54.0; MB1(1,3) = -2.0;
    MB1(2,0) = 36.0; MB1(2,1) = -36.0; MB1(2,2) = -36.0; MB1(2,3) = 36.0;
    MB1(3,0) = -72.0; MB1(3,1) = 216.0; MB1(3,2) = -216.0; MB1(3,3) = 72.0;
  Matrix MB2(3,4); double lambda = sqrt(15.0)/10.0;
    MB2(0,0) = 1.0; MB2(0,1) = lambda; MB2(0,2) = lambda*lambda; MB2(0,3) = lambda*lambda*lambda;
    MB2(1,0) = 1.0; MB2(2,1) = 0.0; MB2(2,2) = 0.0; MB2(2,3) = 0.0;
    MB2(2,0) = 1.0; MB2(2,1) = -lambda; MB2(2,2) = lambda*lambda; MB2(2,3) = -lambda*lambda*lambda;
    Matrix MB(3,4);
    MB = MB2 * MB1/16.0;
  
  Matrix MBT(3,4);   MBT = MB.transpose();
  Matrix Ep(4,2);  // for 4 nodes on an edge
  Matrix p(3,2);   // for Gauss points on an edge
  Matrix EdgeMat(4,4);
  //  DiagMatrix diagRes2Vec(3,3);

  // allocate enough space to create the sparse matrix
  const int MaxContrib = 100*elem.rows() + 16*edges.rows();	
  ColumnVector Si  (MaxContrib);
  ColumnVector Sj  (MaxContrib);
  ColumnVector Sval(MaxContrib);
  
  ColumnVector gVec(nDOF);
  for(int k1 = 0; k1<nDOF; k1++){ gVec(k1) = 0.0; }
  
  /* START THE WORK */	
  int elemTot = elem.rows();
  
  /* loop over all elements */  
  for(int k=0;k<elemTot;k++){ 
    // read element information
    // nodees are called "corners", as I started with the Quad code
    corners(0,0) = nodes((int)elem(k,0)-1,0); //x1
    corners(0,1) = nodes((int)elem(k,0)-1,1); //y1
    corners(1,0) = nodes((int)elem(k,1)-1,0); //x2
    corners(1,1) = nodes((int)elem(k,1)-1,1); //y2
    corners(2,0) = nodes((int)elem(k,2)-1,0); //x3
    corners(2,1) = nodes((int)elem(k,2)-1,1); //y3
    corners(3,0) = nodes((int)elem(k,3)-1,0); //x4
    corners(3,1) = nodes((int)elem(k,3)-1,1); //y4
    corners(4,0) = nodes((int)elem(k,4)-1,0); //x5
    corners(4,1) = nodes((int)elem(k,4)-1,1); //y5
    corners(5,0) = nodes((int)elem(k,5)-1,0); //x6
    corners(5,1) = nodes((int)elem(k,5)-1,1); //y6
    corners(6,0) = nodes((int)elem(k,6)-1,0); //x7
    corners(6,1) = nodes((int)elem(k,6)-1,1); //y7
    corners(7,0) = nodes((int)elem(k,7)-1,0); //x8
    corners(7,1) = nodes((int)elem(k,7)-1,1); //y8
    corners(8,0) = nodes((int)elem(k,8)-1,0); //x9
    corners(8,1) = nodes((int)elem(k,8)-1,1); //y9
    corners(9,0) = nodes((int)elem(k,9)-1,0); //x10
    corners(9,1) = nodes((int)elem(k,9)-1,1); //y10

    // Transformation Matrix - standard triangle to general triangle
    T(0,0) = corners(1,0)-corners(0,0); //x2-x1
    T(0,1) = corners(2,0)-corners(0,0); //x3-x1
    T(1,0) = corners(1,1)-corners(0,1); //y2-y1
    T(1,1) = corners(2,1)-corners(0,1); //y3-y1

    //detT = 2*area;  // or: (x2-x1)(y3-y1)-(x3-x1)(y2-y1)
    detT = T.determinant();  // has to be positive
    
    // fill in diagonal matrices for current element
    for(int gg = 0;gg < 7;gg++){
      diagaxx(gg,gg)  = w(gg)*aV(7*k+gg);
      if (isotropic==1) {
	diagayy(gg,gg) = diagaxx(gg,gg);
	diagaxy(gg,gg) = 0.0;}
      else {
	diagayy(gg,gg) = w(gg)*aV(7*k+gg + nGP);
	diagaxy(gg,gg) = w(gg)*aV(7*k+gg + 2*nGP);
      }
      diagb(gg,gg)  = w(gg)*bV(7*k+gg);
      diagb1(gg,gg) = w(gg)*bxV(7*k+gg);
      diagb2(gg,gg) = w(gg)*byV(7*k+gg);
      wvalVec(gg)   = w(gg)*fV(7*k+gg);
    }

    /* ELEMENT VECTOR */
    /* integration of f*phi */
    elVec = detT*Mtrans*wvalVec;
        
    /* ELEMENT STIFFNESS MATRIX */	
    /* integration of b0*u*phi */	
    elMat = detT*Mtrans*diagb*M;
    /* integration of a*nabla u * nabla phi */
    ttx = T(1,1)*Mxi-T(1,0)*Mnu;   // scalar*7x10 - scalar*7x6 = 7x10
    tty = T(0,0)*Mnu-T(0,1)*Mxi;
    Axx = ttx.transpose()*diagaxx*ttx; // 10x7*7x7*7x7*7x10 = 10x10
    Ayy = tty.transpose()*diagayy*tty;
    Axy = tty.transpose()*diagaxy*ttx;

    elMat += 1/detT*(Axx+Ayy+2.0*Axy); // scalar*(10x10+10x10) = 10x10
    if (convectionFlag){
      /* integration u*b*nabla*phi */ 
      Ab1 = ((+T(1,1)*Mxi.transpose()-T(1,0)*Mnu.transpose())*diagb1*M);
      Ab2 = ((-T(0,1)*Mxi.transpose()+T(0,0)*Mnu.transpose())*diagb2*M);
      elMat -= Ab1+Ab2;
    }

    /* ADD ELEMENT TO GLOBAL STIFFNESS MATRIX */
    for( int iii=0;iii<10;iii++){
      // 0 if not degree of freedom, else the number
      dofs(iii) = (int)n2d((int)elem(k,iii)-1);	
    }
    for (int k1=0; k1<10; k1++){ // loop over nodes
      if (dofs(k1)>0){ // variable value at this node
	for (int k2=0; k2<10; k2++){
	  if (dofs(k2)>0){ // variable value at this node
	    // add values to sparse matrix
	    Si(ptrDOF)   = (int)dofs(k1); 
	    Sj(ptrDOF)   = (int)dofs(k2);
	    Sval(ptrDOF) = elMat(k1,k2);
	    ptrDOF++;
	  }
	  else {// Dirichlet node
	    if (gDscalar){
	      gVec((int)dofs(k1)-1) += elMat(k1,k2)*gDValue;
	    }
	    else {
	      tmpCorner(0) = corners(k2,0);
	      tmpCorner(1) = corners(k2,1);
	      argin(0) = tmpCorner;
	      nargout = 1;
	      res = octave::feval (gDFunc, argin, nargout);
	      gVec((int)dofs(k1)-1) += elMat(k1,k2)*res(0).double_value();
	    }
	  } // endif
	} // endfor k2
	gVec((int)dofs(k1)-1) -= elVec(k1);
      } // endif dofs(k1)>0
    } //endfor k1
  } //endfor k (elements)
  
  /* EDGE CONTRIBUTION */
  //Matrix p(3,2);  // coordinates of the Gauss points
  //Matrix Ep(4,2); // coordinates of the nodes on the edge
  double length;
  ColumnVector vec_diff(2);
  ColumnVector edgeVec(4);
  // insert the edge contributions
  // test if there could be any contribution from edges
  if ((gN1scalar==false)||(gN1Value!=0.0)||(gN2scalar==false)||(gN2Value!=0.0)){
    for(int k=0; k<edges.rows(); k++){ // loop over all edges
      if ((int)edgesT(k)==NEUMANN){  // work with all Neumann edges
	Ep(0,0) = nodes((int)edges(k,0)-1,0);// the four points on the edge
	Ep(0,1) = nodes((int)edges(k,0)-1,1);
	Ep(1,0) = nodes((int)edges(k,1)-1,0);
	Ep(1,1) = nodes((int)edges(k,1)-1,1);
	Ep(2,0) = nodes((int)edges(k,2)-1,0);
	Ep(2,1) = nodes((int)edges(k,2)-1,1);
	Ep(3,0) = nodes((int)edges(k,3)-1,0);
	Ep(3,1) = nodes((int)edges(k,3)-1,1);

	// coodinates of the three Gauss points on the edge
	p(1,0)  = (Ep(0,0)+Ep(3,0))/2.0; // middle point of edge
	p(1,1)  = (Ep(0,1)+Ep(3,1))/2.0;
        vec_diff(0) = sqrt(0.6)*(Ep(3,0)-Ep(0,0))/2.0;
        vec_diff(1) = sqrt(0.6)*(Ep(3,1)-Ep(0,1))/2.0;
	p(0,0) =  p(1,0)+vec_diff(0);
	p(0,1) =  p(1,1)+vec_diff(1);
	p(2,0) =  p(1,0)-vec_diff(0);
	p(2,1) =  p(1,1)-vec_diff(1);

	ColumnVector g(4) ; // Dirichlet values at the four points
	ColumnVector res1Vec(4);
	ColumnVector res2Vec(4);
	length = sqrt(20.0/3.0*(vec_diff(0)*vec_diff(0)+vec_diff(1)*vec_diff(1)));
	// determine the dofs for the four nodes on the boundary
	dofsB(0) = (int)n2d((int)edges(k,0)-1);
	dofsB(1) = (int)n2d((int)edges(k,1)-1);
	dofsB(2) = (int)n2d((int)edges(k,2)-1);
	dofsB(3) = (int)n2d((int)edges(k,3)-1);

	// evaluate the Dirichlet values at the four nodes
	if (gDscalar){	g.fill(gDValue);}
	else{
	  argin(0) = Ep;
	  nargout  = 1;
	  res      = octave::feval(gDFunc, argin, nargout);
	  g        = res(0).column_vector_value();
	}// if (gDscalar)

	// evaluate the values of gN1 at the three Gauss points
	if (gN1scalar){
	  res1Vec.resize(3);
	  res1Vec.fill(gN1Value);
	}
	else{
	  argin(0) = p;
	  nargout  = 1;
	  res      = octave::feval(gN1Func, argin, nargout);
	  res1Vec  = res(0).column_vector_value();
	} // gN1scalar
	
	edgeVec = MBT*(diagwc*res1Vec)*length;

	if (dofsB(0)>0) { gVec((int)dofsB(0)-1) -=  edgeVec(0);}
	if (dofsB(1)>0) { gVec((int)dofsB(1)-1) -=  edgeVec(1);}
	if (dofsB(2)>0) { gVec((int)dofsB(2)-1) -=  edgeVec(2);}
	if (dofsB(3)>0) { gVec((int)dofsB(3)-1) -=  edgeVec(3);}

	// evaluate the values of gN2 at the three Gauss points
	if (gN2scalar){
	  res2Vec.resize(3);
	  res2Vec.fill(gN2Value);
	}
	else{
	  argin(0) = p;
	  nargout  = 1;
	  res      = octave::feval(gN2Func, argin, nargout);
	  res2Vec  = res(0).column_vector_value();
	}// if(gN2scalar)
	
	DiagMatrix diagRes2Vec(res2Vec);
	EdgeMat = length*MBT*(diagwc*diagRes2Vec)*MB;
	
	for (int k1 = 0;k1<4;k1++){
	  if(dofsB(k1)>0){// variable value at this corner, generate equation
	    for(int k2 = 0;k2<4;k2++){
	      if (dofsB(k2)>0){ // variable value at this corner
		// add values to sparse matrix
		Si(ptrDOF)   = (int)dofsB(k1); 
		Sj(ptrDOF)   = (int)dofsB(k2);
		Sval(ptrDOF) = -EdgeMat(k1,k2);
		ptrDOF++;
	      }
	      else {// Dirichlet node at k2
		gVec((int)dofsB(k1)-1) -= EdgeMat(k1,k2)*g(k2);
	      }//if (dofsB(k2)>0
	    }// for k2=0
	  }// dofsB(k1)>0
	}//k1=0
      }// endif ((int)edgesT(k)==NEUMANN){  // work with all Neumann edges
    }//endfor k edges
  }//endif ((gN1scalar==false)||(gN1Value!=0.0)||(gN2scalar==false)||

  /* OUTPUT */
  Si.resize(ptrDOF);  Sj.resize(ptrDOF);  Sval.resize(ptrDOF);
  SparseMatrix sm (Sval,Si,Sj,nDOF,nDOF);
  //  sm.maybe_compress (true);
  retval(0) = sm;
  retval(1) = gVec;
  return retval;
}
