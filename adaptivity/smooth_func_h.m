function refine_flag = smooth_func_h(etpl,etpl_face,refine_flag)
% Smoothing procedure for h-adaptivity.
% Function to ensure new mesh only has 1 hanging node per face.
%
% Input(s):
% etpl        - Element topology structure
% etpl_face   - Element face topology matrix
% refine_flag - Elements to refine in h
%
% Ouput(s):
% refine_flag - Updated Elements to refine in h
%
%  See also SMOOTH_FUNC_POLY.

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% Exit flags
exit_flag=0;
exit_flag_store=zeros(2,1);

% Current level of refinement
refinement_level(:,1)=etpl.tree(:,2)+refine_flag(:,2);

% Interior faces only
etpl_face_i=etpl_face(etpl_face(:,2)>0,1:2);

count=0;
% Infinite loop untill no more changes
while exit_flag==0
    count=count+1;
    % resetting the number of changes
    num_changes=0;
    if count>5000
        fprintf('Error in polynomial smoothing function\n') ;
        pause;
    end
    for j=1:2 % the second loop is the check
        % Looping through all the internal faces
        for i = 1:size(etpl_face_i,1)
            
            % Elements
            el_1 = etpl_face_i(i,1);% element 1
            el_2 = etpl_face_i(i,2);% element 2
            
            % Element polynomial orders
            R_el1=refinement_level(el_1,1);% element 1 polynomial order
            R_el2=refinement_level(el_2,1);% element 2 polynomial order
            
            % Difference in polynomial order
            diff_p=abs(R_el1-R_el2);
            
            % If the difference is greater than 1
            if diff_p>1
                el_change=el_2;
                if R_el1 < R_el2
                    el_change=el_1;
                end
                refinement_level(el_change)=refinement_level(el_change)+1; % increasing the min p element
                num_changes=1; % changes have occured
                %                 fprintf('a change occured\n')
            end
        end
        if num_changes==0
            exit_flag_store(j)=1;
        else
            exit_flag_store=zeros(2,1);
        end
    end
    
    % If two consecutive loops occur with no changes then the mesh is smooth
    if sum(exit_flag_store)==2
        exit_flag=1;
    end
end

refine_flag(:,2)=refinement_level-etpl.tree(:,2);

