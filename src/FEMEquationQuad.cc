// Copyright (C) 2020 Andreas Stahel
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
// Created: 2020-03-30

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/ov-struct.h>

#define NEUMANN   -2
#define DIRICHLET -1

DEFUN_DLD (FEMEquationQuad, args, ,
   "-*- texinfo -*-\n\
@deftypefn {} {} [@var{A},@var{b},@var{n2d}] = FEMEquationQuad(@var{mesh},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2})\n\
\n\
sets up the system of linear equations for a numerical solution of a PDE\n\
using a triangular mesh with elements of order 2\n\
\n\
@verbatim\n\
-div(a*grad u - u*(bx,by))+ b0*u = f         in domain\n\
                               u = gD        on Dirichlet boundary\n\
        n*(a*grad u - u*(bx,by)) = gN1+g2N*u on Neumann boundary\n\
@end verbatim\n\
\n\
parameters:\n\
@itemize\n\
@item @var{mesh} triangular mesh of order 2 describing the domain and the boundary types\n\
@item @var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2}\n\
are the coefficients and functions describing the PDE.\n\
@*Any constant function can be given by its scalar value.\n\
@*The functions @var{a},@var{b0},@var{bx},@var{by} and @var{f} may also be given as vectors\n\
with the values of the function at the Gauss points.\n\
@end itemize\n\
\n\
return values:\n\
@itemize\n\
@item@var{A}, @var{b}: matrix and vector for the linear system to be solved, @var{A}*u-@var{b}=0\n\
@item@var{n2d} is a vectors used to match nodes to degrees of freedom\n\
@* n2d(k)= 0  indicates that node k is a Dirichlet node\n\
@* n2d(k)=nn indicates that the value of the solution at node k is given by u(nn)\n\
@end itemize\n\
@c BEGIN_CUT_TEXINFO\n\
@seealso{FEMEquation, BVP2D, BVP2Dsym, BVP2eig, IBVP2D, I2BVP2D, CreateMeshRect, CreateMeshTriangle}\n\
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
  
  /* INPUT CHECK and Evaluate*/
  string Func;
  //  evaluate function a
  ColumnVector aV; 
  if (args(1).is_string()) {  // function given as string
    Func = args(1).string_value();
    argin(0) = GP;
    nargout  = 1;
    res = octave::feval (Func, argin, nargout);
    aV = res(0).column_vector_value();
  }
  else if(args(1).is_real_scalar()){ //function given as scalar
    aV.resize(GP.rows());
    aV.fill(args(1).double_value());
  }
  else {  // function given by its values
    aV = args(1).column_vector_value();
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
    bV.resize(GP.rows());
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
    bxV.resize(GP.rows());
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
    byV.resize(GP.rows());
    byV.fill(args(4).double_value());
    if (args(4).double_value() !=0){convectionFlag=true;}
  }
  else {  // function given by its values
    bxV = args(3).column_vector_value();
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
    fV.resize(GP.rows());
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
  
  // interpolation matrix for second order elements
  double l1 = (12.0-2.0*sqrt(15.0))/21.0;
  double l2 = (12.0+2.0*sqrt(15.0))/21.0;
  double w1 = (155.0-sqrt(15.0))/2400.0;
  double w2 = (155.0+sqrt(15.0))/2400.0;
  double w3 = 0.1125;

  ColumnVector w(7);  // weights of Gauss points
  w(0)=w1; w(1)=w1; w(2)=w1; w(3)=w2; w(4)=w2; w(5)=w2; w(6)=w3;
  ColumnVector xi(7); // first coordinates of Gauss points
  xi(0)=0.5*l1; xi(1)=1.0-l1; xi(2)=0.5*l1; xi(3)=0.5*l2; xi(4)=1.0-l2; xi(5)=0.5*l2; xi(6)=1.0/3.0;
  ColumnVector nu(7); // second coordinates of Gauss points
  nu(0)=0.5*l1; nu(1)=0.5*l1; nu(2)=1.0-l1; nu(3)=0.5*l2; nu(4)=0.5*l2; nu(5)=1.0-l2; nu(6)=1.0/3.0;
  
  Matrix M(7,6);   // interpolate values of function 
  Matrix Mxi(7,6); // interpolate values of first partial derivative
  Matrix Mnu(7,6); // interpolate values of second partial derivative
  for (int row=0;row<7;row++){
    M(row,0) = (1.0-xi(row)-nu(row))*(1.0-2.0*(xi(row)+nu(row))); 
    M(row,1) = xi(row)*(2*xi(row)-1.0);
    M(row,2) = nu(row)*(2*nu(row)-1.0);
    M(row,3) = 4.0*xi(row)*nu(row);
    M(row,4) = 4.0*nu(row)*(1.0-xi(row)-nu(row));
    M(row,5) = 4.0*xi(row)*(1.0-xi(row)-nu(row));
    
    Mxi(row,0) = -3.0+4.0*(xi(row)+nu(row)); 
    Mxi(row,1) = 4.0*xi(row)-1.0;
    Mxi(row,2) = 0.0;
    Mxi(row,3) = 4.0*nu(row);
    Mxi(row,4) = -4.0*nu(row);
    Mxi(row,5) = 4.0-8.0*xi(row)-4.0*nu(row);
    
    Mnu(row,0) = -3.0+4.0*(xi(row)+nu(row)); 
    Mnu(row,1) = 0.0;
    Mnu(row,2) = 4.0*nu(row)-1.0;
    Mnu(row,3) = 4.0*xi(row);
    Mnu(row,4) = 4.0-4.0*xi(row)-8.0*nu(row);
    Mnu(row,5) = -4.0*xi(row);
  };
  Matrix Mtrans = M.transpose();
  
  
  // diagonal matrices (some will be overwritten later, but have the same size as diagb, therefore they are initialized as follows)
  DiagMatrix diagb(7,7);
  DiagMatrix diaga(7,7);
  DiagMatrix diagb1(7,7);
  DiagMatrix diagb2(7,7);
  // weights for edge contribution	
  ColumnVector wc(3);
  wc(0)=5.0/18.0; wc(1)=8.0/18.0;  wc(2)=5.0/18.0;
  DiagMatrix diagwc(wc);

  // matrices and variables for element stiffness matrix
  Matrix Ab1(6,6);
  Matrix Ab2(6,6);
  Matrix corners(6,2);
  Matrix T(2,2);
  double detT;
  Matrix elMat (6,6);
  ColumnVector elVec (6);
  ColumnVector wvalVec (7);
  Matrix tt(7,6);
  Matrix Ax(6,6);
  Matrix Ay(6,6);
  
  // for assembling of global stiffnes matrix
  RowVector tmpCorner(2); // for jus one corner
  ColumnVector dofs(6);  // degrees of freedom in element
  ColumnVector dofsB(3); // degrees of freedom on edge
  
  // matrices and variables for edge contribution
  double alpha = sqrt(0.6)/2.0;
  Matrix MB(3,3);
  MB(0,0)=alpha*(1.0+2.0*alpha); MB(0,1)=1.0-4.0*alpha*alpha;
                                          MB(0,2)=alpha*(2.0*alpha-1.0);
  MB(1,0) = 0.0;     MB(1,1) = 1.0;               MB(2,1)= 0.0;
  MB(2,0) = MB(0,2); MB(2,1) = MB(0,1);           MB(2,2)= MB(0,0);
  Matrix MBT(3,3);   MBT = MB.transpose();

  Matrix Ep(3,2);  // for 3 points on an edge
  Matrix p(3,2);  // for Gauss points on an edge
  Matrix EdgeMat(3,3);
  //  DiagMatrix diagRes2Vec(3,3);

  // allocate enough space to create the sparse matrix
  const int MaxContrib = 36*elem.rows()+9*edges.rows();	
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

    // Transformation Matrix - standard triangle to general triangle
    T(0,0) = corners(1,0)-corners(0,0); // x2-x1
    T(0,1) = corners(2,0)-corners(0,0); //x3-x1
    T(1,0) = corners(1,1)-corners(0,1); //y2-y1
    T(1,1) = corners(2,1)-corners(0,1); //y3-y1

    //detT = 2*area;  // or: (x2-x1)(y3-y1)-(x3-x1)(y2-y1)
    detT = T.determinant();  // has to be positive
    
    // fill in diagonal matrices for current element
    for(int gg=0;gg<7;gg++){
      diaga(gg,gg)  = w(gg)*aV(7*k+gg);
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
    tt = T(1,1)*Mxi-T(1,0)*Mnu; 	// scalar*7x6 - scalar*7x6 = 7x6
    Ax = tt.transpose()*diaga*tt; // 6x7*7x7*7x7*7x6 = 6x6
    tt = T(0,0)*Mnu-T(0,1)*Mxi;
    Ay = tt.transpose()*diaga*tt;
    elMat += 1/detT*(Ax+Ay); 		// scalar*(6x6+6x6) = 6x6
    if (convectionFlag){
      /* integration u*b*nabla*phi */ 
      Ab1 = ((+T(1,1)*Mxi.transpose()-T(1,0)*Mnu.transpose())*diagb1*M);
      Ab2 = ((-T(0,1)*Mxi.transpose()+T(0,0)*Mnu.transpose())*diagb2*M);
      elMat -= Ab1+Ab2;
    }

    /* ADD ELEMENT TO GLOBAL STIFFNESS MATRIX */
    for( int iii=0;iii<6;iii++){
      // 0 if not degree of freedom, else the number
      dofs(iii) = (int)n2d((int)elem(k,iii)-1);	
    }
    for (int k1=0; k1<6; k1++){ // loop over corners
      if (dofs(k1)>0){ // variable value at this corner
	for (int k2=0; k2<6; k2++){
	  if (dofs(k2)>0){ // variable value at this corner
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
  //  Matrix p(3,2);
  double length;
  ColumnVector vec_diff(2);
  ColumnVector edgeVec(3);
  // insert the edge contributions
  // test if there could be any contribution from edges
  if ((gN1scalar==false)||(gN1Value!=0.0)||(gN2scalar==false)||(gN2Value!=0.0)){
    for(int k=0;k<edges.rows();k++){ // loop over all edges
      if ((int)edgesT(k)==NEUMANN){  // work with all Neumann edges
	Ep(0,0) = nodes((int)edges(k,0)-1,0);// the three points on the edge
	Ep(0,1) = nodes((int)edges(k,0)-1,1);
	Ep(1,0) = nodes((int)edges(k,1)-1,0);
	Ep(1,1) = nodes((int)edges(k,1)-1,1);
	Ep(2,0) = nodes((int)edges(k,2)-1,0);
	Ep(2,1) = nodes((int)edges(k,2)-1,1);
	// coodinates of the three Gauss points on the edge
	p(1,0)  = Ep(1,0); // middle point of edge
	p(1,1)  = Ep(1,1);
        vec_diff(0) = sqrt(0.6)*(Ep(2,0)-Ep(1,0));
        vec_diff(1) = sqrt(0.6)*(Ep(2,1)-Ep(1,1));
	p(0,0) =  p(1,0)-vec_diff(0);
	p(0,1) =  p(1,1)-vec_diff(1);
	p(2,0) =  p(1,0)+vec_diff(0);
	p(2,1) =  p(1,1)+vec_diff(1);

	ColumnVector g(3) ; // Dirichlet values at the three points
	ColumnVector res1Vec(3);
	ColumnVector res2Vec(3);
	length = sqrt(20.0/3.0*(vec_diff(0)*vec_diff(0)+vec_diff(1)*vec_diff(1)));
	// determine the dofs for the three points on the boundary
	dofsB(0) = (int)n2d((int)edges(k,0)-1);
	dofsB(1) = (int)n2d((int)edges(k,1)-1);
	dofsB(2) = (int)n2d((int)edges(k,2)-1);

	// evaluate the Dirichlet values at the three nodes
	if (gDscalar){	g.fill(gDValue);}
	else{
	  argin(0) = Ep;
	  nargout  = 1;
	  res      = octave::feval(gDFunc, argin, nargout);
	  g        = res(0).column_vector_value();
	}// if (gDscalar

	// evaluate the values of gN1 at the three Gauss points
	if (gN1scalar){
	  res1Vec.resize(3);
	  res1Vec.fill(gN1Value);
	}
	else{
	  argin(0) = p;
	  nargout = 1;
	  res     = octave::feval(gN1Func, argin, nargout);
	  res1Vec = res(0).column_vector_value();
	} // gN1scalar
	edgeVec = MBT*(diagwc*res1Vec)*length;
	if (dofsB(0)>0) { gVec((int)dofsB(0)-1) -=  edgeVec(0);}
	if (dofsB(1)>0) { gVec((int)dofsB(1)-1) -=  edgeVec(1);}
	if (dofsB(2)>0) { gVec((int)dofsB(2)-1) -=  edgeVec(2);}

	// evaluate the values of gN2 at the three Gauss points
	if (gN2scalar){
	  res2Vec.resize(3);
	  res2Vec.fill(gN2Value);
	}
	else{
	  argin(0) = p;
	  nargout = 1;
	  res     = octave::feval(gN2Func, argin, nargout);
	  res2Vec = res(0).column_vector_value();
	}// if(gN2sca
	DiagMatrix diagRes2Vec(res2Vec);
	EdgeMat = length*MBT*(diagwc*diagRes2Vec)*MB;
	
	for (int k1=0;k1<3;k1++){
	  if(dofsB(k1)>0){// variable value at this corner, generate equation
	    for(int k2=0;k2<3;k2++){
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
