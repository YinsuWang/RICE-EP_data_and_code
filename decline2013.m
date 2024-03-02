%decline rate of population or TFP calibration
 function [g,delta ] =decline2013(data ,varargin)
 %varargin: min value of g and delta
%   delta+1=g(t-1)/g(t)
%   g(t)+1=data(t)/data(t-1)
g1 = (data(2:size(data,1),:)./data(1:size(data,1)-1,:))-1;
delta1 = (g1(1:size(g1,1)-1,:)./g1(2:size(g1,1),:))-1;
%这里g和dalta都还是矩阵而非向量
delta = mean(delta1,1);
g = mean(g1,1);

%先求均值，再计算g与delta
N = size(data);
for i =1:3
    data1(i,:) = mean(data(i:N-3+i,:),1);
end
%DICE 2013
g1 = (data1(2:size(data1,1),:)./data1(1:size(data1,1)-1,:))-1;
delta1 = (g1(1:size(g1,1)-1,:)./g1(2:size(g1,1),:))-1;
delta = mean(delta1,1);
g = mean(g1,1);

%DICE 2016
% g1 = 1-((data1(2:size(data1,1),:)./data1(1:size(data1,1)-1,:)));
% delta1 = -log((g1(1:size(g1,1)-1,:)./g1(2:size(g1,1),:)));
% delta = mean(delta1,1);
% g = mean(g1,1);

if nargin>=3
    gmin =  varargin{1,1}*ones(1,size(data,2));
    deltamin = varargin{1,2}*ones(1,size(data,2));
    g = max([g;gmin],[],1);
    delta = max([delta;deltamin],[],1);
end