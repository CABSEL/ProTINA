function dfdt=findiff(tpoints,f,tsi)
% clc
% f=[f0,f1,f2,f3,f4];
% tpoints=[0 1 2 3 4];
% tsi=3;

% Finite Difference
% tpoints: numeric vector of time points in increasing order.
%           e.g. tpoints=[0 3 7 12 18 24];
% f: log2(FC) values corresponding to the time points. The number of f
% elements should be equal to the number of time points.
% tsi: the index of the time poinst of interest. The first order
% derivative (df/dt) of the function for f at t=ts will be obatined.

%%
tdiff=tpoints-tpoints(tsi);

deltsind=find(tdiff>0,1,'first');
signdelts=1;
if isempty(deltsind)
    deltsind=find(tdiff<0,1,'last');
    signdelts=-1;
end
delt=tdiff/abs(tdiff(deltsind));
% delt([tsi,deltsind])=[];
Parind=1:1:length(tpoints);
Parind([tsi,deltsind])=[];
numPar=length(Parind);

xmat=[];
ymat=[];
for j=1:1:numPar
    xmat=[xmat;(delt(Parind).^(j+1))];
    ymat=[ymat;(signdelts)^(j+1)];
end
sol=xmat\-ymat;

dfdt=(-(1+sum(sol))*f(:,tsi) + f(:,deltsind) + f(:,Parind)*sol)/((signdelts + delt(Parind)*sol)*abs(tdiff(deltsind)));
