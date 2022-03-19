function [passcheck]=testCheck(simulation,sim_end,testFlag,Er,DG,L2,ndof,nels,passcheck)
% Check if the test problem has been solved correctly.

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/10 17:09:20 $

        testProblemCheck(simulation,sim_end,testFlag,Er,DG,L2,ndof,nels);
        [passcheckNew]=testfileread(DG,L2,Er,ndof,nels,testFlag,simulation,passcheck);
        passcheck=passcheck+passcheckNew;
        if passcheck==8
         fprintf('All tests passed\n');
        end
end
        