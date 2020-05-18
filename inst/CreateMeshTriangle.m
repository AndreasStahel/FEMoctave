## Copyright (C) 2020 Andreas Stahel
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{function file}{}@var{Mesh} = CreateMeshTriangle(@var{name},@var{xy},@var{area},@var{options})
##
##   Generate files with a mesh using the extenal code triangle
##
##parameters:
##@itemize
##@item @var{name} the base filename: the file @var{name}.poly will be generated then triangle will generate files @var{name}.1.* with the mesh
##@item@var{xy} vector containing the coordinates of the nodes forming the outer boundary. The last given node will be connected to the first given node to create a closed curve. Currently no holes can be generated.
##@*The format for @var{xy} is
##       [x1,y1,b1;x2,y2,b2;...;xn,yn,bn] where
##@itemize
##@item xi x-coordinate of node i
##@item yi y-coordinate of node i
##@item bi boundary marker for segment from node i to node i+1
##@itemize
##@item   bi=-1  Dirichlet boundary condition
##@item   bi=-2  Neumann or Robin boundary condition
##@end itemize
##@end itemize
##@item@var{area} the typical area of he individual triangles to be used
##@item@var{options} additional options to be used when calling triangle the options "pa" and the area will be added automatically
##@*Default options are "q" resp. "qpa"
##@*to suppress the verbose information use "Q"
##@end itemize
##
##The information on the mesh generated is written to files and returned in the structure @var{Mesh}, if specified
##@itemize
##@item  The information can then be read and used by
##@*Mesh = ReadMeshTriangle('@var{name}.1');
##@item @var{Mesh} is a a structure with the information about the mesh.
##@*The mesh consists of n_e elements, n_n nodes and n_ed edges.
##@itemize
##@item@var{Mesh.elem} n_e by 3 matrix with the numbers of the nodes forming triangular elements
##@item@var{Mesh.elemArea} n_e vector with the areas of the elements
##@item@var{Mesh.elemT} n_e vector with the type of elements (not used)
##@item@var{Mesh.nodes} n_n by 2 matrix with the coordinates of the nodes
##@item@var{Mesh.nodesT} n_n vector with the type of nodes (not used)
##@item@var{Mesh.edges} n_ed by 2 matrix with the numbers of the nodes forming edges
##@item@var{Mesh.edgesT} n_ed vector with the type of edge
##@item@var{Mesh.GP} n_e*3 by 2 matrix with the coordinates of the Gauss points
##@item@var{Mesh.GPT} n_e*3 vector of integers with the type of Gauss points
##@item@var{Mesh.nDOF} number of DOF, degrees of freedom
##@item@var{Mesh.node2DOF} n_n vector of integer, mapping nodes to DOF
##@end itemize
##@end itemize
##
## Sample call:
##@verbatim
##Mesh = CreateMeshTriangle('Test',[0,-1,-1;1,-1,-2;1,2,-1;0,2,-2],0.01)
##   will create a mesh with 0<=x<=1, -1<=y<=+2
##   and a typical area of 0.01 for each triangle
##Could be read by Mesh = ReadMeshTriangle('Test.1')
##@end verbatim
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitely there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{ReadMeshTriangle, CreateMeshRect, BVP2D, BVP2Dsym, BVP2eig}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2020-03-30

function Mesh = CreateMeshTriangle(name,xy,area,options)

  if ((nargin<3)||(nargin>4)) 
    print_usage()
  endif
  
  opt = '-Q';
  if (nargin==4) opt = options; endif
  
  n = length(xy);

  [fid,msg] = fopen(strcat(name,".poly"),"w");
  fprintf(fid,"# nodes\n");
  fprintf(fid,"%i 2 0 1\n",n);
  fprintf(fid,"%i %16.12f %16.12f %d\n",1,xy(1,1),xy(1,2),max(xy(1,3),xy(n,3)));
  for k = 2:n
    fprintf(fid,"%i %16.15e %16.15e %d\n",k,xy(k,1),xy(k,2),max(xy(k,3),xy(k-1,3)));
  endfor

  fprintf(fid,"# segments\n");
  fprintf(fid,"%i 1\n",n);
  for k = 1:n-1
    fprintf(fid,"%i %i %i %i\n",k,k,k+1,xy(k,3));
  endfor
  fprintf(fid,"%i %i %i %i\n\n",n,n,1,xy(n,3));

  fprintf(fid,"# holes\n");
  fprintf(fid,"0\n",n);

  fclose(fid);
  %% maximal angle 30 degree
  system(['triangle ',opt,sprintf('pq30a%f  ',area),name,'.poly']);
%%  system(['CuthillMcKee -s0 -g ',name,'.1']);

  if nargout >=1
    Mesh = ReadMeshTriangle([name,'.1']);
    eval(['delete ',name,'.1.*'])
  endif
endfunction
