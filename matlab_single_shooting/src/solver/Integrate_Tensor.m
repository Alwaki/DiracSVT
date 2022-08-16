% DESCRIPTION: Uses Runge-Kutta-4 method to integrate dirac equation.

function [ outputFG ] = Integrate_Tensor(xstart, xend, iniF, iniG, k, m, ...
    B, sigmaV0, deltaV0, sigmaR, deltaR, sigmaa, deltaa, tensorV, isospin, Z)

% Declare output array and instanciate
outputFG = [iniF; iniG];

% Step size (TODO: make h argument?)
h=0.001;
step=(xend-xstart)*h;

% RK4 Loop
for i = 1:1:1/h
    x = (i-1)*step+xstart;
    k1 = Tensor_fit(x, outputFG, k, m, B, sigmaV0, deltaV0, sigmaR, deltaR, ...
        sigmaa, deltaa, tensorV, isospin, Z);
    k2 = Tensor_fit(x + step/2, outputFG + k1*step/2, k, m, B, sigmaV0, ...
        deltaV0, sigmaR, deltaR, sigmaa, deltaa, tensorV, isospin, Z);
    k3 = Tensor_fit(x + step/2, outputFG + k2*step/2, k, m, B, sigmaV0, ...
        deltaV0, sigmaR, deltaR, sigmaa, deltaa, tensorV, isospin, Z);
    k4 = Tensor_fit(x + step, outputFG + k3*step, k, m, B, sigmaV0, deltaV0, ...
        sigmaR, deltaR, sigmaa, deltaa, tensorV, isospin, Z);
    outputFG = outputFG + (k1 + 2*k2 + 2*k3 + k4) * step/6;
end


end