% Title:        		DiracSVT
%                
% Authors:        		Alexander Kiessling, Daniel Karlsson, 
%						Yuxin Zhao, Chong Qi
%
% Version:				1.0 (03/2023)	
%
% Project Description:  Numerical solution of the Dirac equation with scalar,
%						vector and tensor potentials
%
% File Description:		Main file of program with user specified settings

%% SETUP

% Cleanup
clc; clear all; close all;

% User specified settings
scenario = 1;
element_state_index = 1;

% Add setup path
addpath("include\utility\");

%% CODE

% Solve for all states
if element_state_index == 0
    for element_state_index = 1:89
        
        % Load data
        run("setup");

        % Call solver routine for given element and scenario
        [B, a0, rvals, FGvals] = dirac_solver(params, B, a0, k_so, Tensor_V);
        
        % Print Energy
        sprintf("%s, B: %f, a0: %f", name, B, a0)
    end

% Solve for a specific state
else
    % Load data
    run("setup");

    % Call solver routine for given element and scenario
    [B, a0, rvals, FGvals] = dirac_solver(params, B, a0, k_so, Tensor_V);
    
    % Print Energy
    sprintf("%s, B: %f, a0: %f", name, B, a0)

    % Plot wavefunction
    plotWF(rvals, FGvals)
end


    

