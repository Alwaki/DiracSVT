% DESCRIPTION: This function returns the boundary conditions when the energy is
% positive. For negative energy, utilize BC.

function [Foutbc, Goutbc, Finbc, Ginbc] = BC_pos(p, B, a0, deltaV0, deltaR, ...
                                            deltaa, sigmaV0, sigmaR, ...
                                            sigmaa, m)
if p.k<0
    Foutbc = -a0*p.xmin^(p.l+2)*(-B + sigmaV0/(1+exp(-sigmaR/sigmaa)))/(p.l+2-p.k);
    Goutbc = a0*p.xmin^(p.l+1);
    Finbc = (SphericalBesselJ(p.l, p.xmax) + SphericalBesselY(p.l, p.xmax));
    Ginbc = sqrt(B^2 + 2*B*m)/(B + 2*m)*(SphericalBesselJ(p.l + 1, p.xmax)...
        + SphericalBesselY(p.l + 1, p.xmax));
else
    Foutbc = a0*p.xmin^p.l;
    Goutbc = a0*p.xmin^(p.l+1)*(2*m+B-deltaV0/(1+exp(-deltaR/deltaa)))/(p.l+p.k+1);
    Finbc = (SphericalBesselJ(p.l, p.xmax) + SphericalBesselY(p.l, p.xmax));
    Ginbc = sqrt(B^2 + 2*B*m)/(B + 2*m)*(SphericalBesselJ(p.l - 1, p.xmax)...
        + SphericalBesselY(p.l - 1, p.xmax));
end

norm = sqrt(hypot(SphericalBesselJ(p.l, p.xmax), SphericalBesselY(p.l, p.xmax)));
Finbc = Finbc/norm;
Ginbc = Ginbc/norm;
