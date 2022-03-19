function el_poly = smooth_func_poly(etpl_face,el_poly)
% Smoothing procedure for p-adaptivity.
% Function to ensure in new mesh polynomial order between
% adjacent elements is only a jump of 1.
%
% Input(s):
% etpl        - Element topology structure
% etpl_face   - Element face topology matrix
% el_poly     - Elements to refine in p
%
% Ouput(s):
% el_poly     - Updated Elements to refine in p
%
%  See also SMOOTH_FUNC_H.

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% Exit flag from infinite loop
exit_flag=0;

% Interior faces only
etpl_face_i=etpl_face(etpl_face(:,2)>0,1:2);
exit_flag_store=zeros(2,1);
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
            p_el1=el_poly(el_1,2);% element 1 polynomial order
            p_el2=el_poly(el_2,2);% element 2 polynomial order
            
            % Difference in polynomial order
            diff_p=abs(p_el1-p_el2);
            
            % If the difference is greater than 1
            if diff_p>1
                el_change=el_2;
                if p_el1 < p_el2
                    el_change=el_1;
                end
                el_poly(el_change,2)=el_poly(el_change,2)+1; % increasing the min p element
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