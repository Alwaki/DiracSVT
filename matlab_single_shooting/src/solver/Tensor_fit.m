% DESCRIPTION: Dirac tensor equation. Adds Coulomb forces 
% due to charge if the nucleon is a proton.

function dFGW = Tensor_fit(x, FG, k, m, B, sigmaV0, deltaV0, ...
    sigmaR, deltaR, sigmaa, deltaa, tensorV, isospin, Z)

% Add charge potential
Col = 0;
if isospin == 1
    if x > sigmaR
        Col = 0.0072923 * Z / x;
    else
        Col = 0.0072923 * Z * (3 * sigmaR^2 - x^2) / (2*sigmaR^3);
    end   
end

% Calculations
sigma = sigmaV0 / (1 + exp((x - sigmaR)/sigmaa)) + Col;      % V+S
delta = deltaV0 / (1 + exp((x - deltaR)/deltaa)) + Col;      % V-S
U = tensorV / (1 + exp((x - sigmaR) / sigmaa));              % WS potential
dFGW = [(-B + sigma) * FG(2) + (k/x - U) * FG(1); (2 * m + B - delta) * FG(1) + (U - k/x) * FG(2)];
end