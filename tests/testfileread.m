function [passcheck]=testfileread(DG,L2,Er,ndof,nels,testFlag,simulation,passcheck)
% Compare computed values against pre-computed values.
% Compares the values from the simulation performed on the 
% users operating system with expected values, read from a test file.
%
% Input(s):
% DG-               - Estimated error in the DG norm
% L2-               - Estimated error in the L2 Norm
% Er-               - Error estimate
% ndof-             - Number of degrees of freedom
% nels-             - Number of elements
% testFlag-         - Flag to indicate which test has been run
% simulation-       - Number of simulations flag
% passcheck-        - Index to mark how many tests are passed (when
%                      passcheck==8, all tests passed)
%
% Output(s):
% passcheck-        - Index to mark how many tests are passed (when
%                      passcheck==8, all tests passed)

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/10 17:09:20 $

DGOpen=fopen('ExpectedDGValues.txt','r');                                   % Open and read expected values of DG norm     
DGArray=fscanf(DGOpen,'%f');

L2Open=fopen('ExpectedL2Values.txt','r');                                   % Open and read expected values of L2 norm 
L2Array=fscanf(L2Open,'%f');

ErOpen=fopen('ExpectedErrorValues.txt','r');                                % Open and read expected values of error estimate 
ErArray=fscanf(ErOpen,'%f');

ndofOpen=fopen('ExpectedNdof.txt','r');                                     % Open and read expected ndof values
ndofArray=fscanf(ndofOpen,'%d');

NelsOpen=fopen('ExpectedNels.txt','r');                                     % Open and read expected nels values
NelsArray=fscanf(NelsOpen,'%d');

tol=1E-10;                                                                  % Tolerance

if testFlag==1 && simulation==1                                             % Test 1, simulation 1 (uniform h adaptivity)
    A(1:3)=DGArray(1:3)-DG';A(4:6)=L2Array(1:3)-L2';A(7:9)=ErArray(1:3)-Er';
    A(10:12)=ndofArray(1:3)-ndof'; A(13:15)=NelsArray(1:3)-nels';           %Subtract calculated values from expected values
    index=zeros(max(size(A)),1);
    
    if all(abs(A(1:3))<tol) && all(abs(A(4:6))<tol) && all(abs(A(7:9))<tol)...
            && all(A(10:12)==0) && all(A(13:15)==0)                         % Test passed if first 9 values in A>tol and ndof and nels are as expected
        fprintf('Test Passed\n')
        passcheck=1;                                                            
    else                                                                    % Test failed else statement
        for i=1:9
            if abs(A(i))>tol                                                % Finding which tests were failed (DG ==> Error)
                index(i)=1;
            end
        end
        
        for i=10:15                                                         % Finding which tests were failed (Ndof ==> Nels)
            if A(i)~=0
                index(i)=1;
            end
        end
    end                                                                     
    
    assert(~any(index(1:3)==1),'Inacceptable results of test 1 with uniform h refinement.Difference between computed values of error in the DG norm and expected values exceed tolerance.\n');
    assert(~any(index(4:6)==1),'Inacceptable results of test 1 with uniform h refinement.Difference between computed values of error in the L2 norm and expected values exceed tolerance.\n');
    assert(~any(index(7:9)==1),'Inacceptable results of test 1 with uniform h refinement.Difference between computed values of the error estimate and expected values exceed tolerance.\n');
    assert(any(index(10:12)==1),'Inacceptable results of test 1 with uniform h refinement.Difference between computed values of ndof and expected values exceed tolerance.\n');
    assert(~any(index(12:15)==1),'Inacceptable results of test 1 with uniform h refinement.Difference between computed values of number of elements and expected values exceed tolerance.\n');
    
    
    %-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
elseif testFlag==1 && simulation==2
    A(1:3)=DGArray(4:6)-DG';A(4:6)=L2Array(4:6)-L2';A(7:9)=ErArray(4:6)-Er';
    A(10:12)=ndofArray(4:6)-ndof'; A(13:15)=NelsArray(4:6)-nels';
    index=zeros(max(size(A)),1);
    
    if all(abs(A(1:3))<tol) && all(abs(A(4:6))<tol) && all(abs(A(7:9))<tol)...
            && all(A(10:12)==0) && all(A(13:15)==0)
        fprintf('Test Passed\n')
        passcheck=1;
    else
        for i=1:9
            if abs(A(i))>tol
                index(i)=1;
            end
        end
        
        for i=10:15
            if A(i)~=0
                index(i)=1;
            end
        end
    end
    
    assert(~any(index(1:3)==1),'Inacceptable results of test 1 with uniform p refinement.Difference between computed values of error in the DG norm and expected values exceed tolerance.\n');
    assert(~any(index(4:6)==1),'Inacceptable results of test 1 with uniform p refinement.Difference between computed values of error in the L2 norm and expected values exceed tolerance.\n');
    assert(~any(index(7:9)==1),'Inacceptable results of test 1 with uniform p refinement.Difference between computed values of the error estimate and expected values exceed tolerance.\n');
    assert(~any(index(10:12)==1),'Inacceptable results of test 1 with uniform p refinement.Difference between computed values of ndof and expected values exceed tolerance.\n');
    assert(~any(index(12:15)==1),'Inacceptable results of test 1 with uniform p refinement.Difference between computed values of number of elements and expected values exceed tolerance.\n');
    
    
    %-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
elseif testFlag==1 && simulation==3
    A(1:3)=DGArray(7:9)-DG';A(4:6)=L2Array(7:9)-L2';A(7:9)=ErArray(7:9)-Er';
    A(10:12)=ndofArray(7:9)-ndof'; A(13:15)=NelsArray(7:9)-nels';
    index=zeros(max(size(A)),1);
        
    if all(abs(A(1:3))<tol) && all(abs(A(4:6))<tol) && all(abs(A(7:9))<tol)...
            && all(A(10:12)==0) && all(A(13:15)==0)
        fprintf('Test Passed\n')
        passcheck=1;
    else
        
        for i=1:9
            if abs(A(i))>tol
                index(i)=1;
            end
        end
        
        for i=10:15
            if A(i)~=0
                index(i)=1;
            end
        end
    end
    
    assert(~any(index(1:3)==1),'Inacceptable results of test 1 with adaptive h refinement.Difference between computed values of error in the DG norm and expected values exceed tolerance.\n');
    assert(~any(index(4:6)==1),'Inacceptable results of test 1 with adaptive h refinement.Difference between computed values of error in the L2 norm and expected values exceed tolerance.\n');
    assert(~any(index(7:9)==1),'Inacceptable results of test 1 with adaptive h refinement.Difference between computed values of the error estimate and expected values exceed tolerance.\n');
    assert(~any(index(10:12)==1),'Inacceptable results of test 1 with adaptive h refinement.Difference between computed values of ndof and expected values exceed tolerance.\n');
    assert(~any(index(12:15)==1),'Inacceptable results of test 1 with adaptive h refinement.Difference between computed values of number of elements and expected values exceed tolerance.\n');
    
    
    %-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif testFlag==1 && simulation==4
    A(1:3)=DGArray(10:12)-DG';A(4:6)=L2Array(10:12)-L2';A(7:9)=ErArray(10:12)-Er';
    A(10:12)=ndofArray(10:12)-ndof'; A(13:15)=NelsArray(10:12)-nels';
    index=zeros(max(size(A)),1);
    if all(abs(A(1:3))<tol) && all(abs(A(4:6))<tol) && all(abs(A(7:9))<tol)...
            && all(A(10:12)==0) && all(A(13:15)==0)
        fprintf('Test Passed\n')
        passcheck=1;
    else
        index=zeros(max(size(A)),1);
        for i=1:9
            if abs(A(i))>tol
                index(i)=1;
            end
        end
        
        for i=10:15
            if A(i)~=0
                index(i)=1;
            end
        end
    end
    
    assert(~any(index(1:3)==1),'Inacceptable results of test 1 with adaptive p refinement.Difference between computed values of error in the DG norm and expected values exceed tolerance.\n');
    assert(~any(index(4:6)==1),'Inacceptable results of test 1 with adaptive p refinement.Difference between computed values of error in the L2 norm and expected values exceed tolerance.\n');
    assert(~any(index(7:9)==1),'Inacceptable results of test 1 with adaptive p refinement.Difference between computed values of the error estimate and expected values exceed tolerance.\n');
    assert(~any(index(10:12)==1),'Inacceptable results of test 1 with adaptive p refinement.Difference between computed values of ndof and expected values exceed tolerance.\n');
    assert(~any(index(12:15)==1),'Inacceptable results of test 1 with adaptive p refinement.Difference between computed values of number of elements and expected values exceed tolerance.\n');
    
    
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif testFlag==1 && simulation==5
    A(1:3)=DGArray(13:15)-DG';A(4:6)=L2Array(13:15)-L2';A(7:9)=ErArray(13:15)-Er';
    A(10:12)=ndofArray(13:15)-ndof'; A(13:15)=NelsArray(13:15)-nels';
    index=zeros(max(size(A)),1);
    if all(abs(A(1:3))<tol) && all(abs(A(4:6))<tol) && all(abs(A(7:9))<tol)...
            && all(A(10:12)==0) && all(A(13:15)==0)
        fprintf('Test Passed\n')
        passcheck=1;
    else
        index=zeros(max(size(A)),1);
        for i=1:9
            if abs(A(i))>tol
                index(i)=1;
            end
        end
        
        for i=10:15
            if A(i)~=0
                index(i)=1;
            end
        end
    end
    
    assert(~any(index(1:3)==1),'Inacceptable results of test 1 with adaptive hp refinement.Difference between computed values of error in the DG norm and expected values exceed tolerance.\n');
    assert(~any(index(4:6)==1),'Inacceptable results of test 1 with adaptive hp refinement.Difference between computed values of error in the L2 norm and expected values exceed tolerance.\n');
    assert(~any(index(7:9)==1),'Inacceptable results of test 1 with adaptive hp refinement.Difference between computed values of the error estimate and expected values exceed tolerance.\n');
    assert(~any(index(10:12)==1),'Inacceptable results of test 1 with adaptive hp refinement.Difference between computed values of ndof and expected values exceed tolerance.\n');
    assert(~any(index(12:15)==1),'Inacceptable results of test 1 with adaptive hp refinement.Difference between computed values of number of elements and expected values exceed tolerance.\n');
    
    
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif testFlag==2
    A(1)=DGArray(16)-DG';A(2)=L2Array(16)-L2';A(3)=ErArray(16)-Er';
    A(4)=ndofArray(16)-ndof'; A(5)=NelsArray(16)-nels';
    index=zeros(max(size(A)),1);
    if abs(A(1))<tol && abs(A(2))<tol && abs(A(3))<tol...
            && A(4)==0 && A(5)==0
        fprintf('Test Passed\n')
        passcheck=1;
    else
        for i=1:3
            if abs(A(i))>tol
                index(i)=1;
            end
        end
        
        for i=4:5
            if A(i)~=0
                index(i)=1;
            end
        end
    end
    
    assert(~index(1)==1,'Inacceptable results of test 2.Difference between computed values of error in the DG norm and expected values exceed tolerance.\n');
    assert(~index(2)==1,'Inacceptable results of test 2.Difference between computed values of error in the L2 norm and expected values exceed tolerance.\n');
    assert(~index(3)==1,'Inacceptable results of test 2.Difference between computed values of the error estimate and expected values exceed tolerance.\n');
    assert(~index(4)==1,'Inacceptable results of test 2.Difference between computed values of ndof and expected values exceed tolerance.\n');
    assert(~index(5)==1,'Inacceptable results of test 2.Difference between computed values of number of elements and expected values exceed tolerance.\n');
    
    
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif testFlag==3
    A=zeros(1,5);
    index=zeros(max(size(A)),1);
    A(1)=DGArray(17)-DG';A(2)=L2Array(17)-L2';A(3)=ErArray(17)-Er';
    A(4)=ndofArray(17)-ndof'; A(5)=NelsArray(17)-nels';
    if abs(A(1))<tol && abs(A(2))<tol && abs(A(3))<tol...
            && A(4)==0 && A(5)==0
        fprintf('Test Passed\n')
        passcheck=1;
    else
        for i=1:3
            if abs(A(i))>tol
                index(i)=1;
            end
        end
        
        for i=4:5
            if A(i)~=0
                index(i)=1;
            end
        end
    end
    
    assert(~index(1)==1,'Inacceptable results of test 3.Difference between computed values of error in the DG norm and expected values exceed tolerance.\n');
    assert(~index(2)==1,'Inacceptable results of test 3.Difference between computed values of error in the L2 norm and expected values exceed tolerance.\n');
    assert(~index(3)==1,'Inacceptable results of test 3.Difference between computed values of the error estimate and expected values exceed tolerance.\n');
    assert(~index(4)==1,'Inacceptable results of test 3.Difference between computed values of ndof and expected values exceed tolerance.\n');
    assert(~index(5)==1,'Inacceptable results of test 3.Difference between computed values of number of elements and expected values exceed tolerance.\n');
    
    
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif testFlag==4
    A=zeros(1,5);
    index=zeros(max(size(A)),1);
    A(1)=DGArray(18)-DG';A(2)=L2Array(18)-L2';A(3)=ErArray(18)-Er';
    A(4)=ndofArray(18)-ndof'; A(5)=NelsArray(18)-nels';
    if abs(A(1))<tol && abs(A(2))<tol && abs(A(3))<tol...
            && A(4)==0 && A(5)==0
        fprintf('Test Passed\n')
        passcheck=1;
    else
        for i=1:3
            if abs(A(i))>tol
                index(i)=1;
            end
        end
        
        for i=4:5
            if A(i)~=0
                index(i)=1;
            end
        end
    end
    
    assert(~index(1)==1,'Inacceptable results of test 4.Difference between computed values of error in the DG norm and expected values exceed tolerance.\n');
    assert(~index(2)==1,'Inacceptable results of test 4.Difference between computed values of error in the L2 norm and expected values exceed tolerance.\n');
    assert(~index(3)==1,'Inacceptable results of test 4.Difference between computed values of the error estimate and expected values exceed tolerance.\n');
    assert(~index(4)==1,'Inacceptable results of test 4.Difference between computed values of ndof and expected values exceed tolerance.\n');
    assert(~index(5)==1,'Inacceptable results of test 4.Difference between computed values of number of elements and expected values exceed tolerance.\n');
    
end

end




