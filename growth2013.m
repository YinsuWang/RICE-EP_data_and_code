
function [data,g] = growth2013(Initial,g0,delta,T)
g=zeros(T+1,size(g0,2));
data = zeros(T+1,size(Initial,2));
g(1,:) = g0;
data(1,:) = Initial;
%DICE 2013
for i =2:T+1
    g(i,:) = g(i-1,:)./(1+delta);
    data(i,:) = data(i-1,:).*(1+g(i,:));
end
%DICE 2016
% for i = 2:T+1
%     g(i,:) = g0.*exp(-delta*(i-1));
%     data(i,:) = data(i-1,:)./(1-g(i-1,:));
% end
end