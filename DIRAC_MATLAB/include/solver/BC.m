% DESCRIPTION: This function returns the boundary conditions when the energy is
% negative. For positive energy, utilize BC_pos.

function [Foutbc, Goutbc, Finbc, Ginbc] = BC(p, B, a0, deltaV0, deltaR, ...
                                            deltaa, sigmaV0, sigmaR, ...
                                            sigmaa, m)
if p.k<0
    Foutbc = -a0*p.xmin^(p.l+2)*(-B+sigmaV0/(1+exp(-sigmaR/sigmaa)))/(p.l+2-p.k);
    Goutbc =  a0*p.xmin^(p.l+1);
else
    Foutbc =  a0*p.xmin^p.l;
    Goutbc =  a0*p.xmin^(p.l+1)*(2*m+B-deltaV0/(1+exp(-deltaR/deltaa)))/(p.l+p.k+1);
end 
    miu=sqrt(-2*m*B-B^2);
    Finbc=-sqrt(-(B)/(2*m+B))*exp(-miu*p.xmax);
    Ginbc=exp(-miu*p.xmax);
end