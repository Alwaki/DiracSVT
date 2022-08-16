%% SINGLE SHOOTING FOR NUMERICAL SOLUTION OF DIRAC'S EQUATION
%
%  Author:          Alexander Wallén Kiessling (akie@kth.se)
%
%  Date:            July 2022
%
%  Description:     This program is intended to calculate the binding
%                   energy of nucleons by solving Dirac's equation. This is
%                   performed with the single shooting method.
%
%  Dependencies:    


%% SETUP

% Cleanup
clc; clear all; close all;

% User specified settings
scenario = 1;
element_state_index = 1;

% Setup acts to set paths, load data, declare variables and enact user settings.
addpath("src\utility\");
run("setup");

%% CODE

% Call solver routine for given element and scenario
[B, a0] = dirac_solver(params, B, a0, k_so, Tensor_V);

    

