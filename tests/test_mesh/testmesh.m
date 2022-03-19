function [coord,etpl,etpl_face]=testmesh(testFlag)
% Generate the mesh for the test problems.
%
% Input(s):
% testFlag     -Flag to determine which problem is running
%
% Ouput(s):
% coord        -Element coordinates
% etpl         -Element topology structure
% etpl_face    -Element face topology

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/10 17:09:20 $

Coord=fopen('Coord Mesh.txt','r');
EtplMat=fopen('EtplMat.txt','r');
EtplPoly=fopen('EtplPoly.txt','r');
EtplTree=fopen('EtplTree.txt','r');
EtplFace1=fopen('EtplFace1.txt','r');
EtplFace2=fopen('EtplFace2.txt','r');
EtplFace3=fopen('EtplFace3.txt','r');
EtplFace4=fopen('EtplFace4.txt','r');

coord=fscanf(Coord,'%f %f',[2 Inf]);
coord=coord';

etplMat=fscanf(EtplMat,'%f %f %f',[3 Inf]);
etpl.mat=etplMat';

etplPoly=fscanf(EtplPoly,'%f %f',[2 Inf]);
etpl.poly=etplPoly';

etplTree=fscanf(EtplTree,'%f %f %f',[3 Inf]);
etpl.tree=etplTree';

if testFlag==1
etplFace=fscanf(EtplFace1,'%f %f %f %f %f %f %f',[7 Inf]);
etpl_face=etplFace';

elseif testFlag==2
etplFace=fscanf(EtplFace2,'%f %f %f %f %f %f %f',[7 Inf]);
etpl_face=etplFace';

elseif testFlag==3
etplFace=fscanf(EtplFace3,'%f %f %f %f %f %f %f',[7 Inf]);
etpl_face=etplFace';

elseif testFlag==4
etplFace=fscanf(EtplFace4,'%f %f %f %f %f %f %f',[7 Inf]);
etpl_face=etplFace';
end
end