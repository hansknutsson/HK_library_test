function [h,im]=gopimage(x,y,z,s)
    
% GOPIMAGE(Z) displays a colour image of the complex matrix Z.
%   The argument is encoded by the gop color table and the
%   magnitude by the intensity. The magnitude is scaled so that the
%   maximum magnitude corresponds to intensity one.
%   If the global variable IMAGEGAMMA exists in the workspace, the
%   image is gamma corrected by that value.
%
% GOPIMAGE(Z), where Z is an MxNx2 real array displays a vector
%   field analogously.
%
% GOPIMAGE(X,Y) does the same thing, but here the components are
%   split into two matrices.
%
% GOPIMAGE(Z,M) and GOPIMAGE(X,Y,M), where M is a scalar, scales
%   the magnitude so that the value M gives maximum intensity.
%
% GOPIMAGE(...,'invert') blends the color determined by the
%   argument with white instead of black, giving images that are
%   better suited for slides.
%
% H=GOPIMAGE(...) returns a handle to the image object.
%
% [H,IM]=GOPIMAGE(...) returns a handle to the image object and the
%   image itself. 
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden

global IMAGEGAMMA	% Global gamma value, defined in workspace
  
% Check the input parameters.
error(nargchk(1,4,nargin))

% Should the image be inverted?
invert=0;
nin=nargin;
switch nin
 case 2
  if ischar(y)
    invert=1;
  end
 case 3
  if ischar(z)
    invert=1;
  end
 case 4
  if ischar(s)
    invert=1;
  end
end
if invert
  nin=nin-1;
end
    
% Is a max level provided?
M=-1;
switch nin
 case 2
  if prod(size(y))==1
    M=y;
    nin=1;
  end
 case 3
  if prod(size(z))==1
    M=z;
    nin=2;
  else
    error('The third parameter must be a scalar.')
  end
end
M=1.0/M; % Convert to normalization factor

if nin==1
  if ndims(x)==2
    u=real(x);
    v=imag(x);
  elseif ndims(x)==3
    if size(x,3)~=2
      error('Z must have a third dimension of size 2.')
    else
      u=x(:,:,1);
      v=x(:,:,2);
    end
  else
    error('Too many dimensions.')
  end
else % nin==2
  if size(x)~=size(y)
    error('X and Y does not have the same size')
  elseif ndims(x)>2
    error('Too many dimensions.')
  else
    u=x;
    v=y;
  end
end

% Get gamma value
if exist('IMAGEGAMMA','var') & isa(IMAGEGAMMA,'double') & ~ ...
      isempty(IMAGEGAMMA) 
  g = IMAGEGAMMA(1);
else
  g = 1.0;
end

w=u+i*v;
gtab=goptab; % Load gop color table
gtab=[gtab;gtab(1,:)]; % Make it cyclic to allow interpolation for all angles
tabangles=2*pi*(0:256)/256;

magw=abs(w);
if (M < 0)
  M = 1.0/max(max(magw)); % Autonormalization
end
magw=min(1,M*magw).^g; % Gamma correct luminance
argw=pi+angle(-w);
argim=reshape(interp1(tabangles,gtab,argw(:)),[size(argw) 3]);
if invert
  gopim=argim.*repmat(magw,[1 1 3])+repmat(1-magw,[1 1 3]);
else
  gopim=argim.*repmat(magw,[1 1 3]);
end
hh=image(gopim);

if nargout>0
  h=hh;
end

if nargout>1
  im=gopim;
end
