function [h] = plot_3D_distribution(coords,varargin)
% --------
% function [h] = plot_3D_distribution(coords,varargin)
% --------
% @ Hans Knutsson   January 2020
%
% Input arguments and default values
% coords                          % A 3D [N,3] vector distribution must be given
% colval      = radi';            % values used for coloring of position markers
% ashell      = [0.1 1];          % alfa - opacity of sphere shells
% rshell      = [eps 1];          % radi of transparent shells
% rhalf       = rad1+0.02;        % radius obove which only positions having x>0 are shown
% ballsize    = 15/ntot^(1/6);    % size of charge marker balls
% shellexp    = [2 2];            % superquadric exponent for transparent shells
% relcolscale = 1;                % 1 => stretch colorspace over min-max values
% shinyballs  = 1;                % 1 => make 3D rendering of each ball
% arrows      = 0;                % 1 => make 3D arrow rendering for each ball
% linewid     = 5;                % width of lines from origin to markers
% linecol     = 0;                % 0 => no lines from origin to markers
% facecol     = 0.12*[1 0.7 1];   % color of 'inner' shell surface
% rotang      = 0;                % rotation angle around z-axis
% metric      = [1 1 1 1 1];      % metric used in simulation [r_wgt a_wgt r_exp r_exp2 a_exp]
% chdens      = [1 1 0 0];        % charge density used in simulation
% rscale      = 1;                % scales ball size
% fignr       = 100;

addpath Util           % Contains various help functions
addpath X_Distribs     % Contains a number of optimized vector distributions

if nargin<1              fprintf('\n A 3D [3,N] vector distribution must be given \n\n');  return; end 

ntot  = size(coords,1);
radi  = sqrt(coords(:,1).^2 + coords(:,2).^2 +coords(:,3).^2);
sradi = sort(radi);
rad1  = sradi(min(8,ntot));   % try to catch first 'shell' with more than 3 orientations

% Set default values for input arguments
colval      = radi';            % values used for coloring of position markers
ashell      = [0.5 0.3];          % alfa - opacity of sphere shells
rshell      = [eps 1];          % radi of transparent shells
rhalf       = 1e6;%rad1+0.02;        % radius obove which only positions having x>0 are shown
ballsize    = 15/ntot^(1/6);    % size of charge marker balls
shellexp    = [2 2];            % superquadric exponent for transparent shells
relcolscale = 1;                % 1 => stretch colorspace over min-max values
shinyballs  = 1;                % 1 => make 3D rendering of each ball
arrows      = 0;                % 1 => make 3D arrow rendering for each ball
linewid     = 5;                % width of lines from origin to markers
linecol     = 0;                % 0 => no lines from origin to markers
facecol     = 0.12*[1 0.7 1];   % color of 'inner' shell surface
rotang      = 0;                % rotation angle around z-axis
metric      = [1 1 1 1 1];      % metric used in simulation [r_wgt a_wgt r_exp r_exp2 a_exp]
chdens      = [1 1 0 0];        % charge density used in simulation
rscale      = 1;                % scales ball size
    
for nvin = 1:2:numel(varargin)
    switch varargin{nvin}
      case 'colval';        colval        = varargin{nvin+1};
      case 'ashell';        ashell        = varargin{nvin+1};
      case 'rshell';        rshell        = varargin{nvin+1};
      case 'rhalf' ;        rhalf         = varargin{nvin+1};
      case 'ballsize';      ballsize      = varargin{nvin+1};
      case 'shellexp';      shellexp      = varargin{nvin+1};
      case 'relcolscale';   relcolscale   = varargin{nvin+1};
      case 'shinyballs';    shinyballs    = varargin{nvin+1};
      case 'arrows';        arrows        = varargin{nvin+1};
      case 'linewidth';     linewid       = varargin{nvin+1};
      case 'linecol';       linecol       = varargin{nvin+1};
      case 'facecol';       facecol       = varargin{nvin+1};
      case 'rotang';        rotang        = varargin{nvin+1};
      case 'metric';        metric        = varargin{nvin+1};
      case 'chdens';        chdens        = varargin{nvin+1};
      case 'rscale';        rscale        = varargin{nvin+1};
      case 'saveCode';      saveCode      = varargin{nvin+1};
      case 'fignr';                         figure(varargin{nvin+1})
      otherwise
        help plot_3D_distribution;
        error([' - Unexpected Option: ' varargin{nvin}])
    end
end
if rhalf<0;     rhalf = rad1+0.02; end 

r_wgt   = metric(1);     % w_r     in MICCAI paper
a_wgt   = metric(2);     % w_phi   in MICCAI paper
r_exp   = metric(3);     % alpha   in MICCAI paper
r_exp2  = metric(4);     % beta    in MICCAI paper
a_exp   = metric(5);     % gamma   in MICCAI paper
chrexp0 = chdens(1);     % r exponent for charge density (r below rmin)
chrexp  = chdens(2);     % r exponent for charge density (r above rmin)
chrmin  = chdens(3);     % rmin in charge density formula
if numel(chdens) < 4
    chdens(4) = 0;       % for backward compatability (due to added cchexp param)   
end
cchexp  = chdens(4);     % cos window exponent (introduced July 9 2018)

R2chd0 = radi.^2.*radi.^chrexp0.*(radi.^2 + chrmin.^2).^(-(chrexp+chrexp0)/2);  % R^2 * Charge density. Must be the same as in main program!
R2chd0 = R2chd0.*cos(radi*pi/2).^cchexp; % add cos window
chnorm = mean(R2chd0);

set(gcf,'Renderer','opengl')
edgeCol    = [1 1 1];    
%faceCol2   = [0.3 0.4 0.5];
faceCol2   = 0.99*[1 1 1];
[x,y,z]    = sphere(80);                % shell resolution
shradi     = (x.^2 + y.^2 + z.^2).^0.5;
colval     = sign(colval).*abs(colval).^(eps+(r_exp + r_exp2)/2); % use the radial metric to modify the color of the balls
bsz0       = 0.008*ballsize;
Nballs     = size(coords,1);
Nsphere    = 24;                        % ball resolution
Ncylinder  = 24;
[xb,yb,zb] = sphere(Nsphere);
%[xc,yc,zc] = cylinder(0.025*[1 1],Ncylinder);
Nshaft     = 30;
[xc,yc,zc] = cylinder(0.0035*linewid*[ones(1,Nshaft) 2.4 1.8 1.2 0.6 0],Ncylinder);
Narrow     = size(xc,1);
Vc         = [xc(:)'; yc(:)'; zc(:)'];


for rang = rotang
    % plot inner shell
    suradi = ((abs(x)+eps).^shellexp(1) + (abs(y)+eps).^shellexp(1) + (abs(z)+eps).^shellexp(1)).^(1/shellexp(1));
    rratio = shradi./suradi;
    x1 = rratio.*rshell(1).*x;
    y1 = rratio.*rshell(1).*y;
    z1 = rratio.*rshell(1).*z;

    mesh(x1,y1,z1,'FaceColor',facecol,'EdgeColor',edgeCol,'Ambientstrength',0.8,'Diffusestrength',0.3,'Specularstrength',0.15,...
         'Specularexponent',20,'FaceAlpha',ashell(1),'EdgeAlpha',0,'FaceLighting','Phong');
    axis equal
    hold on

    % plot outer shell
    suradi = (abs(x).^shellexp(2) + abs(y).^shellexp(2) + abs(z).^shellexp(2)).^(1/shellexp(2));
    rratio = shradi./suradi;
    x2 = rratio.*rshell(2).*x.*(sign(x)*(rshell(2)>=rhalf) + (rshell(2)<rhalf));
    y2 = rratio.*rshell(2).*y;
    z2 = rratio.*rshell(2).*z;

    h=mesh(x2,y2,z2,'FaceColor',faceCol2,'EdgeColor',edgeCol,'Ambientstrength',0.2,'Diffusestrength',0.25,'Specularstrength',0.1,...
           'Specularexponent',500,'FaceAlpha',ashell(2),'EdgeAlpha',0,'FaceLighting','Phong','BackFaceLighting','unlit');
    %get(h)
    hidden off

    %%% plot ball marker
    if relcolscale
        spherecol = abs(colval)-min(abs(colval(:)))+eps;
        spherecol = spherecol/max(spherecol(:));
        spherecol = Hist_eq_cols(spherecol);

    else
        spherecol = abs(colval);
    end

    spherecol = (spherecol./max([1 spherecol]))';
    spherecol = 0.7999*spherecol + 0.2*spherecol.^2;
    spherecol = [spherecol; spherecol];
    r = max(0, 1 - 1.5*spherecol).^0.65;
    g = max(0, 0.65*(1 - abs(2*spherecol - 1)).^0.75);
    b = max(0, min(1, 2*spherecol - 1)).^1.25;

    negcol = ([colval colval] < 0)';
    if 0
        % set negcol balls to black
        r          = (1-negcol).*r;
        g          = (1-negcol).*g;
        b          = (1-negcol).*b;
    else
        % set negcol balls to white
        r          = negcol +(1-negcol).*r;
        g          = negcol +(1-negcol).*g;
        b          = negcol +(1-negcol).*b;
    end
    
    spherecols = [r g b];
    showsphere = (abs(colval) > 0);
    
    linecoln = (linecol./max([1 linecol]))';
    linecolm = 0.7999*linecoln + 0.2*linecoln.^2;
    linecol2 = [linecolm; linecolm];
    r = max(0, 1 - 1.5*linecol2);
    g = max(0, 0.6*(1 - abs(2*linecol2 - 1)).^0.75);
    b = max(0, 2*linecol2 - 1).^1.25;
    linecols = [r g b];
    
    for n = 1:Nballs
       radn = radi(n);
        x0   = cos(rang)*coords(n,1) + sin(rang)*coords(n,2);
        y0   = cos(rang)*coords(n,2) - sin(rang)*coords(n,1);
        z0   = coords(n,3);
        
        if (x0 >= 0) | (radn < rhalf);
            chd = radn.^chrexp0.*(radn.^2 + chrmin.^2).^(-(chrexp+chrexp0)/2);  % Charge density. Must be the same as in main program!
            chd = chd.*max(0.05,cos(radn*pi/2)).^cchexp;                        % add cos window   (0.05 to limit ball size)
            chd = chd/chnorm;                                                   % normalize total ball charge
            bsz = bsz0./max(0.1,chd^(1/3));                                     % ball size inversly proportional to chdr=chargedensity^(1/Dim)

            Xn0   = [x0; y0; z0];
            Xhat0 = Xn0/sqrt(sum(Xn0.^2)+eps);
            XXr   = Xhat0*Xhat0';
            XXa   = eye(3) - XXr; 
            XB    = [xb(:) yb(:) zb(:)]';
            
            % Compute local metric according to metric parameters
            rgrad = sqrt(r_wgt)*r_exp*radn^(r_exp-1);      % must correspond to metric in main program
            agrad = sqrt(a_wgt*a_exp)*radn^(r_exp2-1);     % must correspond to metric in main program
            
            XB  = reshape((XXa/(agrad+eps) + XXr/(rgrad+eps))*XB,3,Nsphere+1,Nsphere+1);
            XBr = sqrt(3*mean(XB(:).^2));
           
            if radn < 1e-4   % take care of ball(s) for r close to 0
               bsz1    = rad1/(XBr+eps); %0.5*chrmin/XBr;
               facealf = 0.25;
            else
                bsz1    = bsz;
                facealf = 0.8;
            end
            bsz2 = bsz1;
  
            if bsz1 > 0.15/(XBr+eps); 
                facealf = 0.5*max(0, 1 - 2*bsz1*XBr).^2;  % to avoid large balls obscuring the sceen
            end
            xr = (bsz1*squeeze(XB(1,:,:)) + x0)*rscale;
            yr = (bsz1*squeeze(XB(2,:,:)) + y0)*rscale;
            zr = (bsz1*squeeze(XB(3,:,:)) + z0)*rscale;
             
            if ballsize > 0
                if shinyballs;
                    h=mesh(xr,yr,zr,'FaceColor',spherecols(n,:),'FaceAlpha',facealf*showsphere(n),'EdgeAlpha',0,'Ambientstrength',0.7,'DiffuseStrength',0.4,...
                           'SpecularStrength',0.6,'Specularexponent',20,'FaceLighting','Phong');

                    if bsz1 ~= bsz2
                        xr2 = bsz2*squeeze(XB(1,:,:)) + x0;
                        yr2 = bsz2*squeeze(XB(2,:,:)) + y0;
                        zr2 = bsz2*squeeze(XB(3,:,:)) + coords(n,3);
                        h=mesh(xr2,yr2,zr2,'FaceColor',spherecols(n,:),'FaceAlpha',0.08,'EdgeAlpha',0,'Ambientstrength',0.5,'DiffuseStrength',0.2,...
                               'SpecularStrength',0.99,'Specularexponent',200,'FaceLighting','Phong');
                    end
                    %colorbar                                               % Does not work
                    %caxis([min(spherecols(n,:)),max(spherecols(n,:))])
                else
                    h= plot3(x0,y0,z0,'ok','markerfacecolor',spherecols(n,:),'MarkerSize',ballsize);
                end
            end 
            lexp = 1 + 1/max(0.1,Nballs^0.25 - 16^0.25);
            if sum(linecol(:)) > 0 && linewid > 0;
                if arrows
                    % Compute rotation matrix
                    r0     = sqrt(x0^2 + y0^2 + z0^2);
                    xh     = x0/r0;
                    yh     = y0/r0;
                    %zh     = z0/r0;
                    vcross = [-yh, xh, 0];  % cross product with cylinder axis ([0 0 1])
                    sina   = sqrt(yh^2 + xh^2);
                    cosa   = sqrt(1 - sina^2)*sign(z0);
                    rotax  = vcross/sina;
                    rx     = rotax(1);
                    ry     = rotax(2);
                    rz     = rotax(3);
                    cosai  = 1 - cosa;
                    R      = cosai*[ rx^2     ,  rx*ry    ,  rx*rz     ;...
                                     rx*ry    ,  ry^2     ,  ry*rz     ;...
                                     rx*rz    ,  ry*rz    ,  rz^2     ]          +...
                                   [ cosa     , -sina*rz  ,  sina*ry   ;...
                                     sina*rz  ,  cosa     , -sina*rx   ;...
                                    -sina*ry  ,  sina*rx  ,  cosa     ];
                    % rotate cylinder
                    RVc = R*Vc;
                    Rxc = reshape(RVc(1,:),[Narrow,Ncylinder+1]);
                    Ryc = reshape(RVc(2,:),[Narrow,Ncylinder+1]);
                    Rzc = reshape(RVc(3,:),[Narrow,Ncylinder+1]);
                    
                    surf(Rxc,Ryc,Rzc,'EdgeAlpha',0,'FaceColor',spherecols(n,:),'Ambientstrength',0.6,'DiffuseStrength',0.3,...
                         'SpecularStrength',0.25,'Specularexponent',2, ...
                         'FaceLighting','Phong')
                else
                    line([0 x0],[0 y0],[0 z0],'Color',linecols(n,:),'Linewidth',linewid);
                end
            end
        end
        hold on
    end

    light('position', [-0.25 0 1])
    light('position', [0.5 0 2],'color',0.7*[1 1 1])
    light('position', [-1 0 -0.5])
    light('position', [1 0 0],'color',0.7*[1 1 1])
    camlight(10,10)
    ma = 1 + bsz2;
    axis(1.02*[-ma ma -ma ma -ma ma])
    axis vis3d
    axis off
    view(0,0)
    %axis tight
    grid on
    set(gca,'cameraviewangle',6.5)
    hold off
    if numel(rotang) > 1
        drawnow
    end
end
