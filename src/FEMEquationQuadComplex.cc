// Copyright (C) 2024 Andreas Stahel
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
// Created: 2025-02-20

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/ov-struct.h>

#define NEUMANN   -2
#define DIRICHLET -1

DEFUN_DLD (FEMEquationQuadComplex, args, ,
   "-*- texinfo -*-\n\
@deftypefn {} {} [@var{A},@var{b}] = FEMEquationQuad(@var{mesh},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2})\n\
\n\
sets up the system of linear equations for a numerical solution of a PDE\n\
using a triangular mesh with elements of order 2.\n\
Complex coefficients are permitted.\n\
\n				     \
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
@item The coefficient @var{a} can als be a symmetric matrix @var{a=[axx,axy;axy,ayy]} given by the row vector @var{[axx,ayy,axy]}.\n\
@* It can be given as row vector or as string with the function name or as nx3 matrix with the values at the Gauss points.\
@end itemize\n					\
\n\
return values:\n\
@itemize\n\
@item@var{A}, @var{b}: matrix and vector for the linear system to be solved, @var{A}*u-@var{b}=0\n\
@end itemize\n\
@c BEGIN_CUT_TEXINFO\n\
@seealso{FEMEquation, FEMEquationQuad, FEMEquationCubic, FEMEquationComplex, FEMEquationCubicComplex, BVP2D, BVP2Dsym, BVP2eig, IBVP2D, I2BVP2D, CreateMeshRect, CreateMeshTriangle}\n\
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
  ComplexColumnVector aV, aVtmp;
  int ii;
  if (args(1).is_string()) {  // function given as string
    Func = args(1).string_value();
    argin(0) = GP;
    nargout  = 1;
    res = octave::feval (Func, argin, nargout);
    //    octave_stdout <<"size of input " << res(0).rows() <<" "<< res(0).columns()<<"\n";
    aV = res(0).complex_column_vector_value();
    if (res(0).columns()==3)
      { isotropic = false;} // store all coefficients in one vector aV
                            // first a11, then a22, then a12
  }
  else if(args(1).length()==1){ //function given as scalar
    aV.resize(nGP);
    aV.fill(args(1).complex_value());
  }
  else {// function given by its values
    //    octave_stdout <<"size of input " << args(1).rows() <<" "<< args(1).columns()<<"\n";
    if(args(1).is_complex_scalar()){ //function given as scalar
      aV.resize(nGP);
      aV.fill(args(1).complex_value());
    }
    else {
      //      octave_stdout <<"not a scalar\n";
      if (args(1).columns()==3)  // anisotropic
      { isotropic = false; // store all coefficients in one vector aV
                           // first a11, then a22, then a12
	if (args(1).rows()==1){// only three values given
	  //	  octave_stdout <<"three values\n";
	  aVtmp = args(1).complex_column_vector_value();
	  aV.resize(3*nGP);
	  //	  octave_stdout << "size of aV "<< aV.rows()<<" "<<aV.columns()<<"\n" ;
	  for(ii = 0; ii < nGP; ii = ii+1){
	    aV(ii)       = aVtmp(0);
	    aV(ii+nGP)   = aVtmp(1);
	    aV(ii+2*nGP) = aVtmp(2);
	  }
	  //octave_stdout << "aVtmp = "<< aVtmp(0)<<" "<<aVtmp(1)<<" "<<aVtmp(2)<<"\n" ;
	  //octave_stdout << "a = "<< aV(1)<<" "<<aV(1+nGP)<<" "<<aV(1+2*nGP)<<"\n" ;

	}
	else {//  all a values in three columns
	  //	  octave_stdout <<"three columns\n";
	  aV = args(1).complex_column_vector_value();
	  //octave_stdout <<"size "<< aV.rows()<<" "<< aV.columns()<<"\n";
	}
      }
      else {//  all a values in three columns
	//	octave_stdout <<"one column\n";
	aV = args(1).complex_column_vector_value();
	//octave_stdout <<"size "<< aV.rows()<<" "<< aV.columns()<<"\n";
      }
    }
  }
  //octave_stdout <<"aV = "<< aV<<"\n";
  //octave_stdout <<"isotropic = "<< isotropic<<"\n";

// evaluate function b
ComplexColumnVector bV; 
if (args(2).is_string()) {  // function given as string
  Func = args(2).string_value();
  nargout  = 1;
  argin(0) = GP;
  res = octave::feval (Func, argin, nargout);
  bV = res(0).complex_column_vector_value();
 }
 else if(args(2).length()==1){ //function given as scalar
   bV.resize(nGP);
   bV.fill(args(2).complex_value());
 }
 else {  // function given by its values
   bV = args(2).complex_column_vector_value();
 }
//cout<<"bV = "<<bV<<"\n";
//cout<<"BLA = "<< args(2).length()<<"\n";

bool convectionFlag = false;
// evaluate function bx
ComplexColumnVector bxV; 
if (args(3).is_string()) {  // function given as string
  convectionFlag = true;
  Func    = args(3).string_value();
  nargout = 1;
  argin(0) = GP;
  res = octave::feval (Func, argin, nargout);
  bxV=res(0).complex_column_vector_value();
 }
 else if(args(3).length()==1){ //function given as scalar
   bxV.resize(nGP);
   bxV.fill(args(3).complex_value());
   if (abs(args(3).complex_value())>0.0){convectionFlag=true;}
 }
 else {  // function given by its values
   bxV = args(3).complex_column_vector_value();
   convectionFlag = true;
 }

// evaluate function by
ComplexColumnVector byV; 
if (args(4).is_string()) {  // function given as string
  Func = args(4).string_value();
  convectionFlag = true;
  nargout  = 1;
  argin(0) = GP;
  res = octave::feval (Func, argin, nargout);
  byV = res(0).complex_column_vector_value();
 }
 else if(args(4).length()==1){ //function given as scalar
   byV.resize(nGP);
   byV.fill(args(4).complex_value());
   if (abs(args(4).complex_value()) >0.0){convectionFlag=true;}
 }
 else {  // function given by its values
   byV = args(4).complex_column_vector_value();
   convectionFlag = true;
 }

//  evaluate function f
ComplexColumnVector fV;
if (args(5).is_string()) {  // function given as string
  Func = args(5).string_value();
  nargout  = 1;
  argin(0) = GP;
  res = octave::feval (Func, argin, nargout);
  fV = res(0).complex_column_vector_value();
 }
 else if(args(5).length()==1){ //function given as scalar
   fV.resize(nGP);
   fV.fill(args(5).complex_value());
 }
 else {  // function given by its values
   fV = args(5).complex_column_vector_value();
 }
 //cout<<"fV = "<<fV<<"\n";
 
//  string gDFunc  = args(6).string_value();
string gDFunc;
Complex gDValue(0.0,0.0);;
bool gDscalar  = true;
if (args(6).is_string()) {  // function given as string
  gDFunc   = args(6).string_value();
  gDscalar = false;
 }
 else {gDValue = args(6).complex_value();
   //cout << "gDValue = "<< gDValue<<"\n";
 }

//  string gN1Func  = args(7).string_value();
string gN1Func;
Complex gN1Value(0.0,0.0);
bool gN1scalar  = true;
if (args(7).is_string()) {  // function given as string
  gN1Func   = args(7).string_value();
  gN1scalar = false;
 }
 else {gN1Value = args(7).complex_value();}

//  string gN2Func  = args(8).string_value();
string gN2Func;
Complex gN2Value(0.0,0.0);
bool gN2scalar  = true;
if (args(8).is_string()) {  // function given as string
  gN2Func   = args(8).string_value();
  gN2scalar = false;
 }
 else {gN2Value = args(8).complex_value();}

/* END INIT */

/* MEMORY ALLOC AND CONST */ 
int ptrDOF = 0;

// interpolation matrix for second order elements
double l1  = (12.0-2.0*sqrt(15.0))/21.0;
double l2  = (12.0+2.0*sqrt(15.0))/21.0;
Complex w1 = (155.0-sqrt(15.0))/2400.0;
Complex w2 = (155.0+sqrt(15.0))/2400.0;
Complex w3 = 0.1125;

ComplexColumnVector w(7);  // weights of Gauss points
w(0)=w1; w(1)=w1; w(2)=w1; w(3)=w2; w(4)=w2; w(5)=w2; w(6)=w3;
ColumnVector xi(7); // first coordinates of Gauss points
xi(0)=0.5*l1; xi(1)=1.0-l1; xi(2)=0.5*l1; xi(3)=0.5*l2; xi(4)=1.0-l2; xi(5)=0.5*l2; xi(6)=1.0/3.0;
ColumnVector nu(7); // second coordinates of Gauss points
nu(0)=0.5*l1; nu(1)=0.5*l1; nu(2)=1.0-l1; nu(3)=0.5*l2; nu(4)=0.5*l2; nu(5)=1.0-l2; nu(6)=1.0/3.0;

ComplexMatrix M(7,6);   // interpolate values of function 
ComplexMatrix Mxi(7,6); // interpolate values of first partial derivative
ComplexMatrix Mnu(7,6); // interpolate values of second partial derivative
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
ComplexMatrix Mtrans = M.transpose();
//cout<<"Mxi = "<<Mxi<<"\n";
//cout<<"Mnu = "<<Mnu<<"\n";
  
  // diagonal matrices (some will be overwritten later, but have the same size as diagb, therefore they are initialized as follows)
ComplexDiagMatrix diagb(7,7);
ComplexDiagMatrix diagaxx(7,7);
ComplexDiagMatrix diagayy(7,7);
ComplexDiagMatrix diagaxy(7,7);
ComplexDiagMatrix diagb1(7,7);
ComplexDiagMatrix diagb2(7,7);
// weights for edge contribution	
ComplexColumnVector wc(3);
//wc(0)=(5.0/18.0,0.0); wc(1)=(8.0/18.0,0.0);  wc(2)=(5.0/18.0,0.0);
 wc(0)=5.0/18.0 ; wc(1)=8.0/18.0 ;  wc(2)=5.0/18.0;
ComplexDiagMatrix diagwc(wc);

// matrices and variables for element stiffness matrix
ComplexMatrix Ab1(6,6);
ComplexMatrix Ab2(6,6);
Matrix corners(6,2);
Matrix T(2,2);
double detT ;
ComplexMatrix elMat (6,6);
ComplexColumnVector elVec (6);
ComplexColumnVector wvalVec (7);
ComplexMatrix ttx(7,6);
ComplexMatrix tty(7,6);
ComplexMatrix Axx(6,6);
ComplexMatrix Ayy(6,6);
ComplexMatrix Axy(6,6);
  
// for assembling of global stiffness matrix
RowVector tmpCorner(2); // for just one corner
ColumnVector dofs(6);  // degrees of freedom in element
ColumnVector dofsB(3); // degrees of freedom on edge

// matrices and variables for edge contribution
double alpha = sqrt(0.6)/2.0;
Matrix MB(3,3);
MB(0,0) = alpha*(1.0+2.0*alpha); MB(0,1)=1.0-4.0*alpha*alpha;
MB(0,2) = alpha*(2.0*alpha-1.0);
MB(1,0) = 0.0;     MB(1,1) = 1.0;               MB(2,1)= 0.0;
MB(2,0) = MB(0,2); MB(2,1) = MB(0,1);           MB(2,2)= MB(0,0);
Matrix MBT(3,3);   MBT = MB.transpose();

Matrix Ep(3,2);  // for 3 points on an edge
Matrix p(3,2);   // for Gauss points on an edge
ComplexMatrix EdgeMat(3,3);
// ComplexDiagMatrix diagRes2Vec(3,3);

// allocate enough space to create the sparse matrix
const int MaxContrib = 36*elem.rows()+9*edges.rows();	
ColumnVector Si  (MaxContrib);
ColumnVector Sj  (MaxContrib);
ComplexColumnVector Sval(MaxContrib);

ComplexColumnVector gVec(nDOF);
 for(int k1 = 0; k1<nDOF; k1++){ gVec(k1) = (0.0,0.0); }

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
  T(0,1) = corners(2,0)-corners(0,0); // x3-x1
  T(1,0) = corners(1,1)-corners(0,1); // y2-y1
  T(1,1) = corners(2,1)-corners(0,1); // y3-y1
  
  //detT = 2*area;  // or: (x2-x1)(y3-y1)-(x3-x1)(y2-y1)
  detT = T.determinant();  // has to be positive

  // fill in diagonal matrices for current element
  for(int gg=0;gg<7;gg++){
    diagaxx(gg,gg)  = w(gg)*aV(7*k+gg);
    if (isotropic==1) {
      diagayy(gg,gg) = diagaxx(gg,gg);
      diagaxy(gg,gg) = (0.0,0.0);}
    else {
      diagayy(gg,gg) = w(gg)*aV(7*k+gg + nGP);
      diagaxy(gg,gg) = w(gg)*aV(7*k+gg + 2*nGP);
    }
    diagb(gg,gg)  = w(gg)*bV(7*k+gg);
    diagb1(gg,gg) = w(gg)*bxV(7*k+gg);
    diagb2(gg,gg) = w(gg)*byV(7*k+gg);
    wvalVec(gg)   = w(gg)*fV(7*k+gg);
  }
  //cout<<"wvalVec = "<< wvalVec <<"\n";
  /* ELEMENT VECTOR */
  /* integration of f*phi */
  elVec = detT*Mtrans*wvalVec;
  //cout <<"detT = "<<detT<<"\n";
  //cout <<"Mtrans = "<<Mtrans<<"\n";
  //cout <<"elVec = "<<elVec<<"\n";
  //cout <<"diagb = "<<diagb<<"\n";
  /* ELEMENT STIFFNESS MATRIX */	
  /* integration of b0*u*phi */	
  elMat = detT*Mtrans*diagb*M;
  /* integration of a*nabla u * nabla phi */
  ttx  = T(1,1)*Mxi-T(1,0)*Mnu; 	// scalar*7x6 - scalar*7x6 = 7x6
  tty  = T(0,0)*Mnu-T(0,1)*Mxi;
  Axx = ttx.transpose()*diagaxx*ttx; // 6x7*7x7*7x7*7x6 = 6x6
  Ayy = tty.transpose()*diagayy*tty;
  Axy = tty.transpose()*diagaxy*ttx;
  //cout<<"ttx="<<ttx<<"\n";
  //cout<<"diagaxx="<<diagaxx<<"\n";
  //cout<<"Axx="<<Axx<<"\n";
  //cout<<"elMat0="<<elMat<<"\n";
  elMat += 1/detT*(Axx+Ayy+2.0*Axy); // scalar*(6x6+6x6) = 6x6
  //cout<<"elMat1="<<elMat<<"\n";
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
	    //cout <<"gDValue"<<gDValue<<"\n";
	    //cout <<"elMat="<<elMat(k1,k2)<<"\n";
	    //cout <<"gVec = "<< gVec;
	  }
	  else {
	    tmpCorner(0) = corners(k2,0);
	    tmpCorner(1) = corners(k2,1);
	    argin(0) = tmpCorner;
	    nargout = 1;
	    res = octave::feval (gDFunc, argin, nargout);
	    gVec((int)dofs(k1)-1) += elMat(k1,k2)*res(0).complex_value();
	  }
	} // endif
      } // endfor k2
      gVec((int)dofs(k1)-1) -= elVec(k1);
    } // endif dofs(k1)>0
  } //endfor k1
 } //endfor k (elements)

/*
// test for complex comparison
 Complex c00(0.0,0.0);
 Complex c01(0.0,1.5);
 Complex c10(1.5,0.0);
 cout<<"c00 = "<<c00<<", c01 = "<<c01<<", c10 = "<<c10<<"\n";
 cout<<"abs(c00) = "<<abs(c00)<<", abs(c01) = "<<abs(c01)<<", abs(c10) = "<<abs(c10)<<"\n";  
 if (c00!=0.0){cout<<"c00!=0.0\n";};
 if (c10!=0.0){cout<<"c10!=0.0\n";};
 if (c01!=0.0){cout<<"c01!=0.0\n";};
 if (c00!=(0.0,0.0)){cout<<"c00!=(0.0,0.0)\n";};
 if (c10!=(0.0,0.0)){cout<<"c10!=(0.0,0.0)\n";};
 if (c01!=(0.0,0.0)){cout<<"c01!=(0.0,0.0)\n";};
 if (abs(c00)>0.0){cout<<"abs(c00)>0\n";};
 if (abs(c10)>0.0){cout<<"abs(c10)>0, abs(c10) = "<<abs(c10)<<"\n";};
 cout<<"c10 = "<<c10<<" abs(c10) = "<<abs(c10)<<"\n";
 if (abs(c01)>0.0){cout<<"abs(c01)>0\n";};
*/
 
/* EDGE CONTRIBUTION */
// Matrix p(3,2);
double length;
ColumnVector vec_diff(2);
ComplexColumnVector edgeVec(3);
// insert the edge contributions
// test if there could be any contribution from edges
 if ((gN1scalar==false)||(abs(gN1Value)>0.0)||(gN2scalar==false)||(abs(gN2Value)>0.0)){
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
      
      ComplexColumnVector g(3) ; // Dirichlet values at the three points
      ComplexColumnVector res1Vec(3);
      ComplexColumnVector res2Vec(3);
      length = sqrt(20.0/3.0*(vec_diff(0)*vec_diff(0)+vec_diff(1)*vec_diff(1)));
      // determine the dofs for the three points on the boundary
      dofsB(0) = (int)n2d((int)edges(k,0)-1);
      dofsB(1) = (int)n2d((int)edges(k,1)-1);
      dofsB(2) = (int)n2d((int)edges(k,2)-1);
      
      // evaluate the Dirichlet values at the three nodes
      if (gDscalar){g.fill(gDValue);}
      else{
	argin(0) = Ep;
	nargout  = 1;
	res      = octave::feval(gDFunc, argin, nargout);
	g        = res(0).complex_column_vector_value();
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
	res1Vec = res(0).complex_column_vector_value();
      } // gN1scalar
      //cout<<"res1Vec = "<<res1Vec<<"\n";
      //cout<<"diagwc = "<<diagwc<<"\n";
      //cout<<"res1Vec = "<<res1Vec<<"\n";
      edgeVec = length*MBT*(diagwc*res1Vec);
      //cout<<"edgeVec = "<<edgeVec<<"\n";
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
	res2Vec = res(0).complex_column_vector_value();
      }// if(gN2scalar)
      ComplexDiagMatrix diagRes2Vec(res2Vec);
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
SparseComplexMatrix sm (Sval,Si,Sj,nDOF,nDOF);
//  sm.maybe_compress (true);
// cout<<gVec;
retval(0) = sm;
retval(1) = gVec;
return retval;
}
