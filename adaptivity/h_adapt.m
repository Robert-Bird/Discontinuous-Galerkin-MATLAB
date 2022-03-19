function [refine_flag]=h_adapt(Er,etpl,etpl_face,d_max)
% Set the refinement flags.
%
% Input(s):
% Er        - Local error for all elements in the mesh
% etpl      - Element topology structure
% etpl_face - Element face topology matrix
% d_max     - d_1 for hp adaptive scheme
%
% Ouput(s):
% refine_flag - elements to refine in h
%
%  See also p_ADAPT, HP_ADAPT..

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

refine_flag = [[1:size(etpl.mat,1)]',zeros(size(etpl.mat,1),1)];
%refinement flag
index=Er>d_max*max(Er);                                                    % Element index for refinement
refine_flag(index,2)=1;                                                    % Refinement flag set
refine_flag = smooth_func_h(etpl,etpl_face,refine_flag);