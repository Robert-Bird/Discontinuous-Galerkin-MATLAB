function [coord,etpl,etpl_face]=crackmesh(~,~,~)
% Opens and reads the .txt files to create the mesh for example problem 3
%
% Input(s):
%
% Ouput(s):
% coord        -Element coordinates
% etpl         -Element topology structure
% etpl_face    -Element face topology

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/21 17:09:20 $

Coord=fopen('crack_Coord.txt','r');                                            % Open .txt files and read
EtplMat=fopen('crack_EtplMat.txt','r');
EtplPoly=fopen('crack_EtplPoly.txt','r');
EtplTree=fopen('crack_EtplTree.txt','r');
EtplFace1=fopen('crack_EtplFace.txt','r');

coord=fscanf(Coord,'%f %f',[2 Inf]);                                       % Generate variables from values in .txt files
coord=coord';

etplMat=fscanf(EtplMat,'%d %d %d',[3 Inf]);
etpl.mat=etplMat';

etplPoly=fscanf(EtplPoly,'%d %d',[2 Inf]);
etpl.poly=etplPoly';

etplTree=fscanf(EtplTree,'%d %d %d',[3 Inf]);
etpl.tree=etplTree';

etplFace=fscanf(EtplFace1,'%d %d %d %d %f %f %d',[7 Inf]);
etpl_face=etplFace';
end