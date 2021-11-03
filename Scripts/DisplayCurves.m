function f = DisplayCurves(Curves,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if size(Curves,1)<size(Curves,2)
    Curves = Curves';
end

if nargin == 1
    f = figure;
elseif nargin == 2
    figure(n);
end

hold on;

s = min(size(Curves,2));
k = ceil(sqrt(s));
for i = 1:s
    subplot(k,k,i);
    hold on;
    plot(Curves(:,i));
end
end

