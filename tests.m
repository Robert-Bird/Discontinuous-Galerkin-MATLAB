% Run test problems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Thomas Wiltshire
% Date:   15/08/2018
% Description: Function to run the test problems (see documentation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable descriptions:
% testFlag-           Flag to indicate which test is being analysed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; path_add; close all;                                                 % Set up problem: clear the workspace, add all folders and their contents to Matlab search path, close all figures
passcheck = 0;

%% Test 1 -----------------------------------------------------------------
testFlag=2;                                                               % testFlag denotes first test problem to run
[sim_end,loop_end,testFlag,E,v,delta2,delta1,av_bc]=testProblem(testFlag);
passcheck_test1 = main_tests(testFlag,sim_end,loop_end,E,v,delta2,delta1,av_bc);
assert(passcheck_test1==3,'Failed test 1');
passcheck=passcheck+passcheck_test1;

%% Test 2 -----------------------------------------------------------------
testFlag=2;                                                                 % testFlag denotes first test problem to run
[sim_end,loop_end,testFlag,E,v,delta2,delta1,av_bc]=testProblem(testFlag);
passcheck_test2 = main_tests(testFlag,sim_end,loop_end,E,v,delta2,delta1,av_bc);
assert(passcheck_test2==1,'Failed test 2');
passcheck=passcheck+passcheck_test2;

%% Test 3 -----------------------------------------------------------------
testFlag=3;                                                                 % testFlag denotes first test problem to run
[sim_end,loop_end,testFlag,E,v,delta2,delta1,av_bc]=testProblem(testFlag);
passcheck_test3 = main_tests(testFlag,sim_end,loop_end,E,v,delta2,delta1,av_bc);
assert(passcheck_test3==1,'Failed test 3');
passcheck=passcheck+passcheck_test3;

%% Test 4 -----------------------------------------------------------------
testFlag=4;                                                                 % testFlag denotes first test problem to run
[sim_end,loop_end,testFlag,E,v,delta2,delta1,av_bc]=testProblem(testFlag);
passcheck_test4 = main_tests(testFlag,sim_end,loop_end,E,v,delta2,delta1,av_bc);
assert(passcheck_test4==1,'Failed test 4');
passcheck=passcheck+passcheck_test4;

if passcheck==8
    fprintf('All tests passed\n');
end



