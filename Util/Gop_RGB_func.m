function RGB = Gop_RGB_func(cim,val,sat,gamma)
% --------------------------------------------------------------------------------------------------------------------
% function RGB = Gop_RGB_func(cim,val,sat,gamma)
% Create Gop style RGB-image from complex input image
% @ Hans Knutsson,  Havsnäs,  Dec 22  2014
%
% cim        complex image 
% val        magnitude              { Optional           real >= 0, default = abs(cim)                               }
% sat        saturation             { Optional           real >= 0, default = 1.          (Can be an image)          }
% gamma      color value exponent   { Optional           real >= 0, default = 0.25                                   }
% RGB        Output RGB image       { Output RGB values  (Truncated at 0 and scaled to max=1) } size = [size(cim),3] }
%
% To show a gop colored image:
% image(RGB)
%
% To plot a gop colored surface:
% surf(abs(cim),RGB)
% ---------------------------------------------------------------------------------------------------------------------

if nargin < 2 || numel(val)   == 0;    val   = abs(cim)/max(abs(cim(:)));   end
if nargin < 3 || numel(sat)   == 0;    sat   = 1;                           end
if nargin < 4 || numel(gamma) == 0;    gamma = 0.25;                        end

val    = val.^gamma;
cangle = angle(cim);                                                           % Complex angle of input image
ca01   = mod(cangle + pi, 2*pi)/2/pi;                                          % ca01 prop to angle and in range [0 1]
Phi    = 2*pi*(0.2*ca01.^2 + 0.8*ca01 - 0.5) + pi/6;                           % Adjust angle to make Gop colors
R      = 0.95*max(0,1-1.7*abs(sin(Phi/2 + 0.3*pi)).^3.5);                      % Red component
G      = 0.9*max(0,1-1.8*abs(sin(Phi/2 + 0.07*pi)).^4);                        % Green component
B      = max(0,1-2.4*abs(sin(Phi/2 - 0.28*pi)).^3);                            % Blue component
dim    = max(find(size(R) > 1));
RGB    = cat(dim+1,R,G,B);                                                         % Color hue
%RGBsz  = size(RGB)
RGB    = bsxfun(@times, val, bsxfun(@plus, bsxfun(@times,sat,RGB), 1-sat));    % Incorporate value and saturation
RGB    = max(0,RGB/max([1; RGB(:)]));                                          % Truncate to 0 and scale to 1 => 0 < RGB < 1

