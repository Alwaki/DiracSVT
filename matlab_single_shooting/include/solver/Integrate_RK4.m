% DESCRIPTION: Uses Runge-Kutta-4 method to integrate dirac equation.

function [ outputFG, rvals, FGvals ] = Integrate_RK4(xstart, xend, iniF, iniG, k, m, ...
    B, sigmaV0, deltaV0, sigmaR, deltaR, sigmaa, deltaa, tensorV, isospin, Z)

% Declare output array and instanciate
outputFG = [iniF; iniG];

% Step size (TODO: make h argument?)
h=0.001;
step=(xend-xstart)*h;

% Storage (radius and wave function)
rvals = zeros(1, 1/h);
FGvals = zeros(2, 1/h);

% RK4 Loop
for i = 1:1:1/h
    x = (i-1)*step+xstart;
    k1 = point_solve(x, outputFG, k, m, B, sigmaV0, deltaV0, sigmaR, deltaR, ...
        sigmaa, deltaa, tensorV, isospin, Z);
    k2 = point_solve(x + step/2, outputFG + k1*step/2, k, m, B, sigmaV0, ...
        deltaV0, sigmaR, deltaR, sigmaa, deltaa, tensorV, isospin, Z);
    k3 = point_solve(x + step/2, outputFG + k2*step/2, k, m, B, sigmaV0, ...
        deltaV0, sigmaR, deltaR, sigmaa, deltaa, tensorV, isospin, Z);
    k4 = point_solve(x + step, outputFG + k3*step, k, m, B, sigmaV0, deltaV0, ...
        sigmaR, deltaR, sigmaa, deltaa, tensorV, isospin, Z);
    outputFG = outputFG + (k1 + 2*k2 + 2*k3 + k4) * step/6;
    
    % Store for plotting
    rvals(i) = x;
    FGvals(1, i) = outputFG(1);
    FGvals(2, i) = outputFG(2);
end


end