% DESCRIPTION: Calculates the spherical Bessel function of order nu at 
% value x, for the first kind of bessel functions.

function js = SphericalBesselJ(eta, z)

if z == 0
    js = besselj(eta, z);
else
    js = sqrt(pi /(2* z)) * besselj(eta + 0.5, z);
end




