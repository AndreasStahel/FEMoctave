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

//  dmat= DiagMatrix(const ColumnVector& a)

DEFUN_DLD (FEMEquation, args, ,
   "-*- texinfo -*-\n\
@deftypefn {} {} [@var{A},@var{b}] = FEMEquation(@var{mesh},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2})\n\
\n\
sets up the system of linear equations for a numerical solution of a PDE\n\
using a triangular mesh with elements of order 1\n\
\n\
@verbatim\n\
-div(a*grad u - u*(bx,by))+ b0*u = f         in domain\n\
                               u = gD        on Dirichlet boundary\n\
        n*(a*grad u - u*(bx,by)) = gN1+g2N*u on Neumann boundary\n\
@end verbatim\n\
\n\
parameters:\n\
@itemize\n\
@item @var{mesh} triangular mesh of order 1 describing the domain and the boundary types\n\
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
@seealso{FEMEquationQuad, FEMEquationCubic, FEMEquationComplex, FEMEquationQuadComplex, FEMEquationCubicComplex, BVP2D, BVP2Dsym, BVP2eig, IBVP2D, I2BVP2D, CreateMeshRect, CreateMeshTriangle}\n\
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

  octave_map arg0 = args(0).map_value ();
  octave_map::const_iterator p1 = arg0.seek ("nodes");
  const Matrix nodes           =  arg0.contents(p1)(0).matrix_value();
  p1 = arg0.seek ("nodesT");
  const ColumnVector nodesT    =  arg0.contents(p1)(0).column_vector_value();
  p1 = arg0.seek ("edges");
  const Matrix edges  =  arg0.contents(p1)(0).matrix_value();
  p1 = arg0.seek ("edgesT");
  const ColumnVector edgesT    =  arg0.contents(p1)(0).column_vector_value();
  p1 = arg0.seek ("elem");
  const Matrix elem  =  arg0.contents(p1)(0).matrix_value();
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
    bxV = res(0).column_vector_value();
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
    convectionFlag=true;
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
  else {gDValue=args(6).double_value();}

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

  int ptrDOF = 0;

  // interpolation matrix for first order elements 
  Matrix M1(3,3);
  M1(0,0) = 2.0/3.0; M1(0,1) = 1.0/6.0; M1(0,2) = 1.0/6.0;
  M1(1,0) = 1.0/6.0; M1(1,1) = 2.0/3.0; M1(1,2) = 1.0/6.0;
  M1(2,0) = 1.0/6.0; M1(2,1) = 1.0/6.0; M1(2,2) = 2.0/3.0;
  
  // allocate enough space to create the sparse matrix
  
  const int MaxContrib = 9*elem.rows()+4*edges.rows();
  ColumnVector Si  (MaxContrib);
  ColumnVector Sj  (MaxContrib);
  ColumnVector Sval(MaxContrib);

  ColumnVector gVec(nDOF);

  for(int k1 = 0; k1<nDOF; k1++){  gVec(k1) = 0.0;}

  Matrix corners(3,2);
  Matrix tmpCorners(1,2);
  double area;

  Matrix M(2,3);
  Matrix elMat (3,3);
  Matrix A(2,2);
  ColumnVector elVec (3);
  ColumnVector dofs(3);

  for(int k=0; k<elem.rows(); k++){ // loop over all elements
    int k3 = 3*k;
    // read element information
    corners(0,0) = nodes((int)elem(k,0)-1,0); 
    corners(0,1) = nodes((int)elem(k,0)-1,1);
    corners(1,0) = nodes((int)elem(k,1)-1,0);
    corners(1,1) = nodes((int)elem(k,1)-1,1);
    corners(2,0) = nodes((int)elem(k,2)-1,0);
    corners(2,1) = nodes((int)elem(k,2)-1,1);

    area = elemArea(k);

    M(0,0) = +corners(2,1)-corners(1,1);
    M(0,1) = +corners(0,1)-corners(2,1);
    M(0,2) = +corners(1,1)-corners(0,1);
    M(1,0) = -corners(2,0)+corners(1,0);
    M(1,1) = -corners(0,0)+corners(2,0);
    M(1,2) = -corners(1,0)+corners(0,0);

    Matrix bMat(3,3);
    for(int ii=0;ii<=2;ii++){
      for(int jj=0;jj<=2;jj++){
	bMat(ii,jj) = M1(ii,jj)*bV(k3+ii)*area/3.0;
      }
    }

    if (isotropic ==1){
      // factor* M'*M
      double factor = (aV(k3)+aV(k3+1)+aV(k3+2))/(12.0*area);
      elMat = factor * M.transpose()*M;
    }
    else{ //anisotropic case
      A(0,0) = (aV(k3)+aV(k3+1)+aV(k3+2));
      A(1,1) = (aV(k3+nGP)+aV(k3+1+nGP)+aV(k3+2+nGP));
      A(1,0) = (aV(k3+2*nGP)+aV(k3+1+2*nGP)+aV(k3+2+2*nGP));
      A(0,1) = A(1,0);
      elMat = 1.0/(12.0*area)*M.transpose()*A*M;
    }

    elMat += M1*bMat;
    elVec(0) = fV(k3)  *area/3.0;
    elVec(1) = fV(k3+1)*area/3.0;
    elVec(2) = fV(k3+2)*area/3.0;
    elVec    = M1*elVec;

    Matrix convMat(3,3);
    double SignArea;
    if (convectionFlag){ // compute only if convection term is there
      convMat(0,0) = -(bxV(k3)  *M(0,0) + byV(k3)  *M(1,0));
      convMat(1,0) = -(bxV(k3)  *M(0,1) + byV(k3)  *M(1,1));
      convMat(2,0) = -(bxV(k3)  *M(0,2) + byV(k3)  *M(1,2));
      convMat(0,1) = -(bxV(k3+1)*M(0,0) + byV(k3+1)*M(1,0));
      convMat(1,1) = -(bxV(k3+1)*M(0,1) + byV(k3+1)*M(1,1));
      convMat(2,1) = -(bxV(k3+1)*M(0,2) + byV(k3+1)*M(1,2));
      convMat(0,2) = -(bxV(k3+2)*M(0,0) + byV(k3+2)*M(1,0));
      convMat(1,2) = -(bxV(k3+2)*M(0,1) + byV(k3+2)*M(1,1));
      convMat(2,2) = -(bxV(k3+2)*M(0,2) + byV(k3+2)*M(1,2));
      SignArea = ((corners(1,0)-corners(0,0))*(corners(2,1)-corners(0,1))-
	          (corners(1,1)-corners(0,1))*(corners(2,0)-corners(0,0)));
      if (SignArea>0) SignArea = 1.0/6.0;
      else SignArea = -1.0/6.0;

      elMat -= SignArea*convMat*M1;
    }
    // now element matrix and vector are generated
    // insert elements in  global matrix and vector
    dofs(0) = (int)n2d((int)elem(k,0)-1);
    dofs(1) = (int)n2d((int)elem(k,1)-1);
    dofs(2) = (int)n2d((int)elem(k,2)-1);


    for (int k1=0;k1<3;k1++){ // loop over corners
      if (dofs(k1)>0){ // variable value at this corner
	for (int k2=0;k2<3;k2++){
	  if (dofs(k2)>0){ // variable value at this corner
	    // add values to sparse matrix
	    Si(ptrDOF)   = (int)dofs(k1); 
	    Sj(ptrDOF)   = (int)dofs(k2);
	    Sval(ptrDOF) = elMat(k1,k2);
	    ptrDOF++;}
	  else {// Dirichlet node
	    if (gDscalar){
	      gVec((int)dofs(k1)-1) += elMat(k1,k2)*gDValue;
	    }
	    else {
	      //Matrix tmpCorners(1,2);
	      tmpCorners(0,0) = corners(k2,0);
	      tmpCorners(0,1) = corners(k2,1);
	      argin(0) = tmpCorners;
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
  
  // insert the edge contributions
  // test if there could be any contribution from edges
  if ((gN1scalar==false)||(gN1Value!=0.0)||(gN2scalar==false)||(gN2Value!=0.0)){
    double factp = (1.0+1.0/sqrt(3.0))/2.0;
    double factn = (1.0-1.0/sqrt(3.0))/2.0;
    for(int k = 0; k<edges.rows(); k++){ // loop over all edges
      if ((int)edgesT(k)==NEUMANN){  // work with all Neumann edges
	corners(0,0) = nodes((int)edges(k,0)-1,0); 
	corners(0,1) = nodes((int)edges(k,0)-1,1);
	corners(1,0) = nodes((int)edges(k,1)-1,0); 
	corners(1,1) = nodes((int)edges(k,1)-1,1);
	Matrix points(2,2);
	points(0,0) = factp*corners(0,0)+factn*corners(1,0);
	points(0,1) = factp*corners(0,1)+factn*corners(1,1);
	points(1,0) = factn*corners(0,0)+factp*corners(1,0);
	points(1,1) = factn*corners(0,1)+factp*corners(1,1);
	argin(0) = points;
	ColumnVector res1Vec;
	ColumnVector res2Vec;
	double alpha  = (1.0-1.0/sqrt(3.0))/2.0;
	double length =
	  sqrt( (corners(1,0)-corners(0,0))*(corners(1,0)-corners(0,0))+
		(corners(1,1)-corners(0,1))*(corners(1,1)-corners(0,1)) )/2.0;
	dofs(0) = (int)n2d((int)edges(k,0)-1);
	dofs(1) = (int)n2d((int)edges(k,1)-1);

	if ((gN1scalar==false)||(gN1Value!=0.0)){
	  // evaluate the first function
	  if (gN1scalar){
	    res1Vec.resize(2);
	    res1Vec.fill(gN1Value);
	  }
	  else{
	    nargout = 1;
	    res     = octave::feval(gN1Func, argin, nargout);
	    res1Vec = res(0).column_vector_value();
	  }
	  ColumnVector edgeVec(2);
	  edgeVec(0) = ((1.0-alpha)*res1Vec(0)+alpha*res1Vec(1)) * length;
	  edgeVec(1) = (alpha*res1Vec(0)+(1.0-alpha)*res1Vec(1)) * length;
	  if (dofs(0)>0) { gVec((int)dofs(0)-1) -=  edgeVec(0);} 
	  if (dofs(1)>0) { gVec((int)dofs(1)-1) -=  edgeVec(1);}
	}
	if ((gN2scalar==false)||(gN2Value!=0.0)){
	  // evaluate the second function
	  if (gN2scalar){
	    res2Vec.resize(2);
	    res2Vec.fill(gN2Value);
	  }
	  else{
	    nargout = 1;
	    res     = octave::feval(gN2Func, argin, nargout);
	    res2Vec = res(0).column_vector_value();
	  }
	  Matrix EdgeMat(2,2);
	  EdgeMat(0,0) = ((1.0-alpha)*(1.0-alpha)*res2Vec(0)
			  +alpha*alpha*res2Vec(1))*length;
	  EdgeMat(0,1) = ((1.0-alpha)*alpha*(res2Vec(0)+res2Vec(1)))*length;
	  EdgeMat(1,0) = EdgeMat(0,1);
	  EdgeMat(1,1) = ((1.0-alpha)*(1.0-alpha)*res2Vec(1)
			  +alpha*alpha*res2Vec(0))*length;
	  for (int k1=0;k1<2;k1++){
	    if(dofs(k1)>0){// variable value at this corner, generate equation
	      for(int k2=0;k2<2;k2++){
		if (dofs(k2)>0){ // Neumann node, variable value at this corner
		  // add values to sparse matrix
		  Si(ptrDOF)   = (int)dofs(k1); 
		  Sj(ptrDOF)   = (int)dofs(k2);
		  Sval(ptrDOF) = -EdgeMat(k1,k2);
		  ptrDOF++;
		}
		else {// Dirichlet node at k2
		  if (gDscalar){
		    gVec((int)dofs(k1)-1)  -= EdgeMat(k1,k2)*gDValue;
		  }
		  else {
		    Matrix tmpCorner(1,2);
		    tmpCorner(0,0) = corners(k2,0);
		    tmpCorner(0,1) = corners(k2,1);
		    argin(0) = tmpCorner;
		    nargout  = 1;
		    res      = octave::feval (gDFunc, argin, nargout);
		    gVec((int)dofs(k1)-1) -= EdgeMat(k1,k2)*res(0).double_value();
		  }// if gDscalar
		}//else Dirichlet node
	      }// k2
	    }// dofs(k1)
	  }//k1
	}// if gN2scalar
      }// endif edgesT==0
    }//endfor k edges
  }//endif edge contributions

  Si.resize(ptrDOF);  Sj.resize(ptrDOF);  Sval.resize(ptrDOF);
  SparseMatrix sm (Sval,Si,Sj,nDOF,nDOF);
  
  retval(0) = sm;
  retval(1) = gVec;

  return retval;
}
