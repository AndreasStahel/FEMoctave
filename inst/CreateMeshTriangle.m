## Copyright (C) 2022 Andreas Stahel
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
##   Generate files with a mesh  with linear elements using the external code triangle
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
##@item   bi = -1  Dirichlet boundary condition
##@item   bi = -2  Neumann or Robin boundary condition
##@end itemize
##@end itemize
##@item@var{area} the typical area of he individual triangles to be used
##@item@var{options} additional options to be used when calling triangle the options "pa" and the area will be added automatically
##@*Default options are "q" resp. "qpa"
##@*to suppress the verbose information use "Q"
##@end itemize
##
##The information on the mesh generated is written to files and returned in the structure @var{Mesh}, if the return argument is provided.
##@itemize
##@item  The information can then be read and used by
##@*Mesh = ReadMeshTriangle('@var{name}.1');
##@item @var{Mesh} is a structure with the information about the mesh.
##@*The mesh consists of n_e elements, n_n nodes and n_ed edges.
##@itemize
##@item@var{Mesh.type} a string with the type of triangle: linear, quadratic or cubic
##@item@var{Mesh.elem} n_e by 3 (or 6/10) matrix with the numbers of the nodes forming triangular elements
##@item@var{Mesh.elemArea} n_e vector with the areas of the elements
##@item@var{Mesh.elemT} n_e vector with the type of elements (not used)
##@item@var{Mesh.nodes} n_n by 2 matrix with the coordinates of the nodes
##@item@var{Mesh.nodesT} n_n vector with the type of nodes (not used)
##@item@var{Mesh.edges} n_ed by 2 (or 3/4) matrix with the numbers of the nodes forming edges
##@item@var{Mesh.edgesT} n_ed vector with the type of edge
##@item@var{Mesh.GP} n_e*(3/7) by 2 matrix with the coordinates of the Gauss points
##@item@var{Mesh.GPT} n_e*(3/7) vector of integers with the type of Gauss points
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
## Created: 2022-11-22

function Mesh = CreateMeshTriangle(name,xy,area,varargin)

  if ((nargin<3))
    print_usage()
  endif
  opt = '-Q';
  CuthillMcKee = 0;
  DeleteFiles = 1;
  
  PointCounter = 0; Points = []; Areas = [];
  HoleCounter  = 0; HoleSize = 0;
  HoleLength   = [];  HoleBorder = []; HolePoint = [];
  SegmentCounter = 0; SegmentBorder = [];
  SegmentLength = []; 
  if (~isempty(varargin))
    for cc = 1:length(varargin)
      switch tolower(varargin{cc}.name)
        case {'meshsize'}
	  NewPoints = varargin{cc}.where;
	  Points = [Points; NewPoints];
	  Areas = [Areas;varargin{cc}.area];
	  PointCounter++;
	case {'hole'}
	  NewBorder  = varargin{cc}.border;
	  HoleBorder = [HoleBorder;NewBorder];
	  HolePoint  = [HolePoint;varargin{cc}.point];
	  HoleLength = [HoleLength;size(NewBorder,1)];
	  HoleCounter++;
	case {'segment'}
	  NewBorder = varargin{cc}.border;
	  SegmentBorder = [SegmentBorder;NewBorder];
	  SegmentLength = [SegmentLength;size(NewBorder,1)];
	  SegmentCounter++;
	case {'option'}
	  switch tolower(varargin{cc}.option)
	    case 'cuthillmckee'
	      CuthillMcKee = varargin{cc}.value;
	    case 'deletefiles'
	      DeleteFiles = varargin{cc}.value;
	    otherwise
	      error('Invalid optional argument, %s.%s',varargin{cc}.name,varargin{cc}.option)
	  endswitch
	otherwise
	  error('Invalid optional argument, %s',varargin{cc}.name);
      endswitch % switch
    endfor % for
  endif % if

  n = length(xy);
  HoleSize = size(HoleBorder,1);
  [fid,msg] = fopen(strcat(name,".poly"),"w");
  fprintf(fid,"# nodes\n");
  fprintf(fid,"%i 2 0 1\n",n+HoleSize+size(SegmentBorder,1));
  fprintf(fid,"%i %16.12f %16.12f %d\n",1,xy(1,1),xy(1,2),max(xy(1,3),xy(n,3)));
  for k = 2:n
    fprintf(fid,"%i %16.15e %16.15e %d\n",k,xy(k,1),xy(k,2),max(xy(k,3),xy(k-1,3)));
  endfor % k
  HoleLimits = 1+[0;cumsum(HoleLength)];
  if HoleCounter>0
    for hc = 1:HoleCounter
      fprintf(fid,"%i %16.15e %16.15e %d\n",n+HoleLimits(hc),...
	      HoleBorder(HoleLimits(hc),1),HoleBorder(HoleLimits(hc),2),...
	      max(HoleBorder(HoleLimits(hc),3),HoleBorder(HoleLimits(hc+1)-1,3)));
      for jj = HoleLimits(hc)-1+[2:HoleLength(hc)];
	fprintf(fid,"%i %16.15e %16.15e %d\n",n+jj,HoleBorder(jj,1),
		HoleBorder(jj,2),max(HoleBorder(jj,3),HoleBorder(jj-1,3)));
      endfor %% jj
    endfor %% hc
  endif %% HoleCounter
  SegmentLimits = [0;cumsum(SegmentLength)];
  if SegmentCounter>0
    for sc = 1:SegmentCounter
      for jj = 1:SegmentLength(sc);
	fprintf(fid,"%i %16.15e %16.15e %d\n",n+HoleSize+SegmentLimits(sc)+jj,...
		SegmentBorder(SegmentLimits(sc)+jj,1),SegmentBorder(SegmentLimits(sc)+jj,2),SegmentBorder(SegmentLimits(sc)+jj,3));
      endfor% jj
    endfor % sc
  endif%% SegmentCounter
  fprintf(fid,"# segments\n");
  if size(SegmentLength,1) == 0
    fprintf(fid,"%i 1\n",n+HoleSize);
  else
    fprintf(fid,"%i 1\n",n+HoleSize+sum(SegmentLength)-SegmentCounter);
  endif
  for k = 1:n-1
    fprintf(fid,"%i %i %i %i\n",k,k,k+1,xy(k,3));
  endfor
  fprintf(fid,"%i %i %i %i\n",n,n,1,xy(n,3));

  for hc = 1:HoleCounter
    for jj = 1:HoleLength(hc)-1;
      base = n+HoleLimits(hc)+jj-1;
      fprintf(fid,"%i %i %i %i\n",base,base,base+1,HoleBorder(HoleLimits(hc)+jj,3))
    endfor
    base = n+HoleLimits(hc+1)-1;
    fprintf(fid,"%i %i %i %i\n",base,base,n+HoleLimits(hc),HoleBorder(HoleLimits(hc),3))
  endfor

  if SegmentCounter>0
    cc = 1;
    for sc = 1:SegmentCounter
      for jj = 1:SegmentLength(sc)-1;
	fprintf(fid,"%i %i %i %i\n",n+HoleSize+cc,...
		n+HoleSize+SegmentLimits(sc)+jj,...
		n+HoleSize+SegmentLimits(sc)+jj+1,...
		SegmentBorder(SegmentLimits(sc)+jj,3))
	cc++;
      endfor% jj
    endfor %% sc
  endif
 
  fprintf(fid,"# holes\n");
  fprintf(fid,"%i\n",HoleCounter);
  for hc = 1:HoleCounter
    fprintf(fid,"%i %16.15e %16.15e\n",hc,HolePoint(hc,1), HolePoint(hc,2));
  endfor
  
  if PointCounter>0
      fprintf(fid,"# area markers\n");
    fprintf(fid,"%i\n",PointCounter);
    for jj = 1:PointCounter
      fprintf(fid,"%i %f %f 0 %f \n",jj,Points(jj,1), Points(jj,2),Areas(jj));
    endfor
  endif
  
  %% maximal angle 30 degree
  if PointCounter>0
    command = ['triangle ',opt,sprintf('pq30a '),name,'.poly'];
  else
    command = ['triangle ',opt,sprintf('pq30a%f ',area),name,'.poly'];
  endif

  fprintf(fid,"# generate mesh by : %s\n",command);
  fclose(fid);
  
  system(command);
  if CuthillMcKee
    system(['./CuthillMcKee -s0 ',name,'.1']);
  endif
  
  if nargout >=1
    Mesh = ReadMeshTriangle([name,'.1']);
    if DeleteFiles
      eval(['delete ',name,'.1.*'])
    endif
  endif
  Mesh.type = 'linear';
endfunction
