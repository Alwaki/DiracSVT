% DESCRIPTION: Calculates the spherical Bessel function of order nu at 
% value x, for the second kind of bessel functions.

function js = SphericalBesselY(eta, z)

if z == 0
    js = bessely(eta, z);
else
    js = sqrt(pi /(2* z)) * bessely(eta + 0.5, z);
end




