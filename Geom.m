function [x,y]=Geom(bs,s) % fid=wgeom(g,'Geom');
%GEOM	Gives geometry data for the Geom PDE model.
%
%   NE=GEOM gives the number of boundary segments
%
%   D=GEOM(BS) gives a matrix with one column for each boundary segment
%   specified in BS.
%   Row 1 contains the start parameter value.
%   Row 2 contains the end parameter value.
%   Row 3 contains the number of the left-hand regions.
%   Row 4 contains the number of the right-hand regions.
%
%   [X,Y]=GEOM(BS,S) gives coordinates of boundary points. BS specifies the
%   boundary segments and S the corresponding parameter values. BS may be
%   a scalar.

nbs=8;

if nargin==0,
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 0 0 0 0 0 % start parameter value
  1 1 1 1 1 1 1 1 % end parameter value
  2 2 2 2 2 2 2 2 % left hand region
  1 1 1 1 0 0 0 0 % right hand region
];

bs1=bs(:)';

if find(bs1<1 | bs1>nbs),
  error(message('pde:wgeom:NonExistBoundSeg'))
end

if nargin==1,
  x=d(:,bs1);
  return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 & n==1,
  bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) | n~=size(s,2),
  error(message('pde:wgeom:BsSizeError'));
end

if ~isempty(s),

% boundary segment 1
ii=find(bs==1);
if length(ii)
x(ii)=(0.5-(0.5))*(s(ii)-d(1,1))/(d(2,1)-d(1,1))+(0.5);
y(ii)=(-0.5-(0.5))*(s(ii)-d(1,1))/(d(2,1)-d(1,1))+(0.5);
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=(-0.5-(0.5))*(s(ii)-d(1,2))/(d(2,2)-d(1,2))+(0.5);
y(ii)=(-0.5-(-0.5))*(s(ii)-d(1,2))/(d(2,2)-d(1,2))+(-0.5);
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=(-0.5-(-0.5))*(s(ii)-d(1,3))/(d(2,3)-d(1,3))+(-0.5);
y(ii)=(0.5-(-0.5))*(s(ii)-d(1,3))/(d(2,3)-d(1,3))+(-0.5);
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=(0.5-(-0.5))*(s(ii)-d(1,4))/(d(2,4)-d(1,4))+(-0.5);
y(ii)=(0.5-(0.5))*(s(ii)-d(1,4))/(d(2,4)-d(1,4))+(0.5);
end

% boundary segment 5
ii=find(bs==5);
if length(ii)
x(ii)=1*cos(1.9634954084936207*s(ii)+(-3.6651914291880918))+(0);
y(ii)=1*sin(1.9634954084936207*s(ii)+(-3.6651914291880918))+(0);
end

% boundary segment 6
ii=find(bs==6);
if length(ii)
x(ii)=1*cos(1.439896632895322*s(ii)+(-1.7016960206944711))+(0);
y(ii)=1*sin(1.439896632895322*s(ii)+(-1.7016960206944711))+(0);
end

% boundary segment 7
ii=find(bs==7);
if length(ii)
x(ii)=1*cos(1.439896632895322*s(ii)+(-0.26179938779914913))+(0);
y(ii)=1*sin(1.439896632895322*s(ii)+(-0.26179938779914913))+(0);
end

% boundary segment 8
ii=find(bs==8);
if length(ii)
x(ii)=1*cos(1.4398966328953215*s(ii)+(1.1780972450961729))+(0);
y(ii)=1*sin(1.4398966328953215*s(ii)+(1.1780972450961729))+(0);
end

end
