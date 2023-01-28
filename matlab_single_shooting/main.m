%% SINGLE SHOOTING FOR NUMERICAL SOLUTION OF DIRAC'S EQUATION
%
% 
% Project:        Shell evolution of the dirac equation
%                 
% Authors:        Alexander W. Kiessling
%                 (2022-2023)
% 
% Description:    Main file of program with user specified settings
% 


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
[B, a0, rvals, FGvals] = dirac_solver(params, B, a0, k_so, Tensor_V);

% Print Energy
sprintf("B: %f, a0: %f", B, a0)

% Plot wavefunction
plot(rvals, FGvals(1, :), rvals, FGvals(2,:))

    

