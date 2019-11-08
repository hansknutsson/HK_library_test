function mygopimage(z, maxval ,gamma, scale)
%    
% function mygopimage(Z, maxval, gamma) displays a colour image of the 
% complex matrix Z.
%
% z      -  A 2D array or a wa (kernel generator array)
%
% maxval -  The magnitude is scaled so that maxval gives maximum
%           intensity. If maxval is negative, the image is auto
%           normalized. (optional, default=1).
%
% gamma  -  Sets the global variable IMAGEGAMMA = gamma
%           this value affects successive plots and prints.
%           (optional, default sqrt(0.5));
%
% scale  - Scale the true image size with this factor (Optional,
%          default 1)
%
% Author: Mats Andersson, matsa@imt.liu.se
%

default_maxval = 1;
default_gamma = sqrt(0.5);
default_scale = 1;

switch nargin
    case 1
	maxval = default_maxval;
	 gamma = default_gamma;
	 scale = default_scale;
     case 2
	 gamma = default_gamma;
	 scale = default_scale;
     case 3
     	 scale = default_scale;
     case 4
     otherwise
        help mygopimage
        return
end     

if maxval <= 0 
    maxval = -1;   % compatibility to earlier version
end

if isa(z, 'wa')     
    z = getdata(z);
end

z = squeeze(z);

global IMAGEGAMMA
IMAGEGAMMA = gamma;
shrinkwrap(scale, gopimage(z, maxval))
figtitle(inputname(1))

