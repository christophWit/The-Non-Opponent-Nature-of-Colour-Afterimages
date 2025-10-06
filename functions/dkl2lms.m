%% dkl2lms
function [lms, M] = dkl2lms(dkl, lms0, monlms)
% Inverts lms2dkl by inverting the matrix.
% 2025.04.23

L0 = lms0(1);
M0 = lms0(2);
S0 = lms0(3);

% MATRIX THAT SCALES RELATIVE TO WP (L0, M0, S0)
M = [...
    1/(L0+M0),  1/(L0+M0),      0;...
    1/L0,       -1/M0,      0;...
    -1/(L0+M0), -1/(L0+M0), 1/S0];

if nargin > 2
    [maxL, maxS] = map2mongamut(M, lms0, monlms);
    M(2,:) = M(2,:)/maxL;
    M(3,:) = M(3,:)/maxS;
end

dlms = M \ dkl';
dlms = dlms';

% INCREMENTS TO LMS:
lms = dlms+ones(size(dlms,1),1)*lms0;
lms = round(lms,12);
