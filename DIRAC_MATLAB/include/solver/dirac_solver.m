% DESCRIPTION: Main solver which loops RK4 routine until convergence or
% divergence of solution.

function [B, a0, rvals, FGvals] = dirac_solver(p, B, a0, k_so, tensorV)

Gap = @(x,y) (x-y);    %% define gap function to retain similarity of codes

% Determine mass and potentials
A = p.N + p.Z;
if p.isospin == 1
     m = 938.27208816;
     sigmaV0 = p.V0 * (1 + p.kappa * (p.N - p.Z) / A);
     if p.scenario == 1
        deltaV0 = -p.lambda * sigmaV0; 
     elseif p.scenario == 2
        deltaV0 = -p.lambda * p.V0 * (1 - p.kappa * (p.N-p.Z) / A); 
     elseif p.scenario == 3
        deltaV0 = -p.lambda * p.V0 * (1 - k_so * (p.N-p.Z) / A);
     end
else
    m = 939.56542052;
    sigmaV0 = p.V0 * (1 - p.kappa*(p.N-p.Z)/A);
    if p.scenario == 1
        deltaV0 = -p.lambda * sigmaV0; 
    elseif p.scenario == 2
        deltaV0 = -p.lambda * p.V0 * (1 + p.kappa*(p.N-p.Z)/A); 
    elseif p.scenario == 3
        deltaV0 = -p.lambda * p.V0 * (1 + k_so*(p.N-p.Z)/A);
    end
end

sigmaR = p.r0 * A^(1/3);
deltaR = p.Rls * A^(1/3);
sigmaa = p.a;
deltaa = p.a;
if p.scenario==3
    deltaa = p.als;
end

% Set solver tolerances and step
h = 0.0001;
error = 1;
convergence_threshold = 0.0001;
iterations = 0;

% Enter solver loop
while (error > convergence_threshold) && abs(B) < 100 && iterations < 200

    % Boundary conditions at zero point 
    if B < 0
        [Foutbc, Goutbc, Finbc, Ginbc] = BC(p, B, a0, deltaV0, deltaR, ...
                                            deltaa, sigmaV0, sigmaR, ...
                                            sigmaa, m);
    else
        [Foutbc, Goutbc, Finbc, Ginbc] = BC_pos(p, B, a0, deltaV0, deltaR, ...
                                            deltaa, sigmaV0, sigmaR, ...
                                            sigmaa, m);
    end

    [inFG, rvals, FGvals] = Integrate_RK4(p.xmax, p.xmatch, Finbc, Ginbc, p.k, m, ...
    B, sigmaV0, deltaV0, sigmaR, deltaR, sigmaa, deltaa, tensorV, p.isospin, p.Z);

    [outFG, outrvals, outFGvals] = Integrate_RK4(p.xmin, p.xmatch, Foutbc, Goutbc, p.k, m, ...
    B, sigmaV0, deltaV0, sigmaR, deltaR, sigmaa, deltaa, tensorV, p.isospin, p.Z);

    rvals = [outrvals flip(rvals)];
    FGvals = [outFGvals flip(FGvals,2)];
    
    % Boundary conditions at infinity point
    B1 = B + B*h;
    if B1 < 0
        [Foutbc, Goutbc, Finbc, Ginbc] = BC(p, B1, a0, deltaV0, deltaR, ...
                                            deltaa, sigmaV0, sigmaR, ...
                                            sigmaa, m);
    else
        [Foutbc, Goutbc, Finbc, Ginbc] = BC_pos(p, B1, a0, deltaV0, deltaR, ...
                                            deltaa, sigmaV0, sigmaR, ...
                                            sigmaa, m);
    end

    dBinFG = Integrate_RK4(p.xmax, p.xmatch, Finbc, Ginbc, p.k, m, ...
    B1, sigmaV0, deltaV0, sigmaR, deltaR, sigmaa, deltaa, tensorV, p.isospin, p.Z);

    dBoutFG = Integrate_RK4(p.xmin, p.xmatch, Foutbc, Goutbc, p.k, m, ...
    B1, sigmaV0, deltaV0, sigmaR, deltaR, sigmaa, deltaa, tensorV, p.isospin, p.Z);
    
    dGFdB    = (Gap(dBoutFG, dBinFG) - Gap(outFG, inFG)) / (B*h);
    da0outGF = outFG * (1 + h);
    dGFda0   = (Gap(da0outGF, inFG) - Gap(outFG, inFG)) / (a0*h);

    M = [dGFdB dGFda0];
    Cold = (M*[B;a0])-Gap(outFG,inFG);
    Old = [B;a0];
    New = M\Cold;
    B = New(1);
    a0 = New(2);
    error = (norm(New-Old));

    iterations = iterations + 1;

end

% Remove divergence errors, simply give these high values
if abs(B) > 100 || isnan(B) || ~isreal(B)
    B = 100;
end