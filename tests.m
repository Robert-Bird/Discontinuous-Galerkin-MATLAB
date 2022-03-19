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
global testFlag

%% Test 1
testFlag = 1;
results = runtests('main_tests.m');
disp(table(results))

%% Test 2
testFlag = 2;
results = runtests('main_tests.m');
disp(table(results))

%% Test 3
testFlag = 3;
results = runtests('main_tests.m');
disp(table(results))

%% Test 4
testFlag = 4;
results = runtests('main_tests.m');
disp(table(results))