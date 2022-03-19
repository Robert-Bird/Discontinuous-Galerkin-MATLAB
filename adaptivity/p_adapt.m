function etpl=p_adapt(Er,etpl,d_max,d_min)
% Application of the adaptive p refinement.
% 
% Input(s):
% Er        - Local error
% etpl      - Element topology structure
% d_max     - d_1 for hp adaptive scheme
% d_min     - d_2 for hp adaptive scheme
%
%
% Ouput(s):
% etpl      - Updated element topology structure
%
%  See also H_ADAPT, HP_ADAPT.

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

Er_max=max(Er);                                                            % Max error
index=(((Er_max*d_min)<Er) .* (Er<=(Er_max*d_max)))==1;                   % p-adaptive strategy 3
etpl.poly(index,2)=etpl.poly(index,2)+1;                                   % Changing polynomial order
