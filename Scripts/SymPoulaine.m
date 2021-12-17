function [Poulaine] = SymPoulaine(HalfPoulaine,Side)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    Side = 'R';
else
    error('Input error');
end
Poulaine = [HalfPoulaine, HalfPoulaine];

mid = fix(size(HalfPoulaine,1)/2);

switch Side
    case 'R'
        Poulaine(:,4:6) = [[Poulaine(1:mid,), Poulaine(:,), Poulaine(:,)] ; [Poulaine(:,), Poulaine(:,), Poulaine(:,)]];
    case 'L'
        Poulaine(:,1:3) = [HalfPoulaine, HalfPoulaine];
end

end

