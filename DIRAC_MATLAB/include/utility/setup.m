% DESCRIPTION: This script is intended to act as a cleanup for the main file.
% Any clutter is distributed here. In sequential order of this script, this
% includes user defined settings, folder pathing, loading file data,
% setting data specific parameters, declaring variables as well as
% allocating space for data structures.

% Add folder structure pathing
addpath("include\")
addpath("include\solver")
addpath("include\utility")
addpath("files\")

% Load data from file
data = readtable("data.xlsx", "VariableNamingRule","preserve");
parameters = readtable("parameters.xlsx", "VariableNamingRule","preserve");

% Read scenario parameters from table
V0 = parameters{scenario,2};
kappa = parameters{scenario,3};
lambda = parameters{scenario,4};
r0 = parameters{scenario,5};
a = parameters{scenario,6};
Rls = parameters{scenario,7};
als = parameters{scenario,8};
k_so = parameters{scenario,9};
Tensor_V = parameters{scenario,10};

% Read element parameters from table
isospin = data{element_state_index, 3};
N = data{element_state_index, 4};
Z = data{element_state_index, 5};
l = data{element_state_index, 6};
k = data{element_state_index, 7};
B = data{element_state_index, 8};
xmin = data{element_state_index, 9};
xmax = data{element_state_index, 10};
xmatch = data{element_state_index, 11};
a0 = data{element_state_index, 12};

% Create parameter container
params = ParameterContainer.createObj(V0, kappa, lambda, r0, a, Rls, ...
                                        als, isospin, N, Z, l, k, ...
                                        xmin, xmax, xmatch, scenario);
