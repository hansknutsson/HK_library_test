function [Vwgts, ErrSPH, A, W] = Optimize_SPH_filters(varargin)
% --------
% function [Vwgts, ErrSPH, A, W] = Optimize_SPH_filters(varargin)
% --------
% Finds the amplitude of the projection of a set of 3 dimensional vectors on spherical harmonics (SPH) functions.
% Only a 'single' vector set should be entered. Using double ([x -x]) sets will work but be slower.
%
% @ Hans Knutsson   January 2019
%
% Inputs and default values:
% --------------------------
% Distr:  [3,N] ([1,1])       - Input vector array (or number of precomputed vectors)
%                              (N = 6,10,16,21,30,46,61,102,120,200,300,400,500,600 or 800 exist,
%                               the closest number will be picked)
%
% Optional input arguments (varargin):
%
% Vwgts0   = [];   Cell      - SPH weight vectors. Cell size = (1, number of SPH degrees).
%                              If Vwgts0 is entered these values are used, i.e no optimization is done.
% cert     = 1;    [0,1]     - ones(1,N), 0 => lost sample
% rcosexp  = 2.8;  Real      - Sets cos^exp type weight function factor in optimization.
% wrW      = 2;    Real      - Sets wrW in (SPHdegree + d0)^(-wrW) type weight function factor in optimization.
% d0       = 6;    Real      - Sets d0 in (SPHdegree + d0)^(-wrW) type weight function factor in optimization.
% rcosexpS = 0.7;  Real      - Sets cos^exp type test SPH spectrum factor.
% wrS      = 1;    Real      - Sets wrS in (SPHdegree + 1)^(-wrS) type test SPH spectrum factor.
% c0frac   = 0;    Real      - Approximate fraction of lost measurements (0 to 1)
% ppars    = 0;    [0,1]     - 1 => print current parameters,  0 => no print.
% peval    = 0;    [0,1]     - 1 => print evaluation numbers,  0 => no print.
% plotlev  = 0;    Integer   - 0 => no plots,  1 => a few plots,  2 = many plots,  3 = all plots.
% rendSPH  = 0;    Integer   - 0< => render color coded SPH weights up to order 2*(rendSPH-1)
% fignr1   = 5000; Integer   - Figure numbers start at fignr1.
% saveCode = 0;    [0,1]     - 1 => The code that is beeing executed will be saved.
%
% Outputs:
% --------
% Vwgts:  [1,N] real array  - Optimized weights for each vector

addpath Util           % Contains various help functions
addpath X_Distribs     % Contains a number of optimized vector distributions

% Set default values for input arguments
vecs     = 61;
Vwgts0   = [];
cert     = [];
rcosexp  = 2.8;
wrW      = 2;                                   
d0       = 6;                                   
rcosexpS = 2.8/4;
wrS      = 1;
c0frac   = 0;
ppars    = 1;
peval    = 1;                                   
plotlev  = 0;   
rendSPH  = 0;
relcol   = 0;
Npdeg    = 2;
fignr1   = 5000; 
saveCode = 0;       

optimize = 1;  % will be set to zero if Vwgts are entered
for nvin = 1:2:numel(varargin)
    switch varargin{nvin}
      case 'Distr';       vecs       = varargin{nvin+1};
      case 'Vwgts';       Vwgts0     = varargin{nvin+1};  optimize   = 0;
      case 'cert';        cert       = varargin{nvin+1};
      case 'rcosexp';     rcosexp    = varargin{nvin+1};
      case 'wrW';         wrW        = varargin{nvin+1};
      case 'd0';          d0         = varargin{nvin+1};
      case 'rcosexpS';    rcosexpS   = varargin{nvin+1};
      case 'wrS';         wrS        = varargin{nvin+1};
      case 'c0frac';      c0frac     = varargin{nvin+1};
      case 'ppars';       ppars      = varargin{nvin+1};
      case 'peval';       peval      = varargin{nvin+1};
      case 'plotlev';     plotlev    = varargin{nvin+1};
      case 'rendSPH';     rendSPH    = varargin{nvin+1};
      case 'relcol';      relcol     = varargin{nvin+1};
      case 'Npdeg';       Npdeg      = varargin{nvin+1};
      case 'fignr1';      fignr1     = varargin{nvin+1};
      case 'saveCode';    saveCode   = varargin{nvin+1};
      otherwise
        help Optimize_SPH_filters;
        error([' - Unexpected Option: ' varargin{nvin}])
    end
end

if numel(vecs) == 1;
    Nvec   = vecs;
    Nexist = [6 10 16 21 30 46 61 102 120 200 300 400 500 600 800];
    Nvec1  = Nexist(find(abs(Nexist - Nvec) - min(abs(Nexist - Nvec)) == 0));  % pick closest existing distribution
    if min(abs(Nexist - Nvec)) > 0; fprintf('\n! Picking closest precomputed vector distribution: Nvec = %3i \n',Nvec1); end
    switch Nvec1
      case 6;   dName = 'X_006_Qdistr_03641__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat';
      case 10;  dName = 'X_010_Qdistr_03640__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat';      
      case 16;  dName = 'X_016_Qdistr_03633__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat';
      case 21;  dName = 'X_021_Qdistr_03629__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat'; % (used for ISMRM 2019)
      case 30;  dName = 'X_030_Qdistr_03631__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat';
      case 46;  dName = 'X_46_SHELL.mat';   % (very good)
      case 61;  dName = 'X_061_Qdistr_03632__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat';
      case 102; dName = 'X_102_Qdistr_03634__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat';
      case 120; dName = 'X_120_SHELL.mat';  % (very good)      
      case 200; dName = 'X_200_oneSHELL_2.00_2.00_2.00_1.00_1.00_J0_S1.mat'; % (very good)
      case 300; dName = 'X_300_Qdistr_03635__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat';
      case 400; dName = 'X_400_Qdistr_03655__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat'
      case 500; dName = 'X_500_Qdistr_03645__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat';
      case 600; dName = 'X_600_Qdistr_03652__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat';
      case 800; dName = 'X_800_Qdistr_03662__1.00_2.00_8.00_1.00_0.00_J0_S1_1.00_1.50_0.03_0.00.mat';
      %otherwise; fprintf('\n No precomputed %3i vector distribution found \n',Nvec); return
    end
    load(dName)
    vecs = rivec1(:,1:2:end);      % don't use diametrical pairs
end
Nvec1 = size(vecs,2);        

% Set measurement certainties if not given
if numel(cert) == 0;
    cert      = ones(1,Nvec1);
    Nc0       = randperm(Nvec1);
    nc0       = Nc0(1:round(min(1,c0frac)*Nvec1));  % randomly selected zero certainties
    cert(nc0) = 0;
end

%cert(1,1) = 0;                % loose one
%cert(1,5) = 0;                % loose one good for X_120 
%cert(1,2:end) = 0;            % leave one
%cert(1,4:end) = 0;            % leave three
%cert(1,11:end) = 0;           % leave ten
%cert(1,[1:100]) = 0;          % loose 100
%cert(1,[1:50]) = 0;           % loose 50
%cert(1,[27:56]) = 0;          % good for X_300
%cert(1,1:20) = 0;             % good for X_200
%cert(1,26:37) = 0;            % good for X_120  (used for ISMRM 2019)
%cert(1,[25:44]) = 0;          % good for X_102
%cert(1,[5  12:16]) = 0;       % good for X_61
%cert(1,[1 2 5 6]) = 0;        % good for X_46
%cert(1,[1 2 5:10]) = 0;       % good for X_46
%cert(1,[21 22]) = 0;          % good for X_46
%cert(1,[3 4 7 8 13 14]) = 0;  % good for X_30
%cert(1,[3 4 7]) = 0;          % good for X_30
%cert(1,[3 4]) = 0;            % good for X_21
%cert(1,[2 5]) = 0;            % good for X_21  (used for ISMRM 2019)
%cert(1,:) = abs(abs(x(3,:))-0.6) <0.25;

maxdeg  = 2*round(0.75*(sqrt(1+3*Nvec1*8/9)-1));  % aprox 3 X above SPH 'Nyqusist'  (must be above SPH Nyq)
Nlost   = sum(cert <= 0);
Nvecok  = Nvec1 - Nlost;

useVwgts0 = numel(Vwgts0) > 0;

% just to test rotation invariance of measures
rp = 0;
if rp ~= 0
    x = [cos(rp) 0 -sin(rp);
         0       1     0;
         sin(rp) 0 cos(rp)]*x;
end

xrad    = sqrt(sum(vecs.^2,1));
xhat    = bsxfun(@rdivide, vecs, xrad);
xhat    = bsxfun(@times, xhat, xrad > 1e-4); % smaller radius considered = 0
degrees = 0:2:maxdeg;
ndegree = [0:maxdeg/2];
degref  = 1.5*(sqrt(1+(Nvecok-1)*8/9)-1);            % Nvecok SPH 'Nyqusist'
maxdegp = min(round(0.75*maxdeg)+4,maxdeg);%2*ceil(0.75*(sqrt(1+1.5*Nvec1*8/9)-1)) + 2;   % more than 1.5 X above SPH 'Nyqusist' 

cosexp0 = ((maxdeg+eps)/(degref+eps))^2 - 1;
%cosexp0old = (maxdeg + 1)*(maxdeg + 2)/2/sqrt(Nvec1*Nvecok)
cosexp  = rcosexp*cosexp0;
cos1    = cos(pi/2*ndegree/(ndegree(end) + 1));
coswgt  = cos1.^cosexp;
dnum    = 2*degrees + 1;
dnum0   = 2*d0 + 1;
degwgt  = (dnum0./(dnum0 + dnum - 1)).^(wrW/2);
sphwgt  = coswgt.*degwgt;                   % Optimization weighting function in SPH domain

cosexpS = rcosexpS*cosexp0;
coswgtS = cos1.^cosexpS;
degwgtS = (1./dnum).^(wrS/2); 
spectS  = coswgtS.*degwgtS;                 % Test signal spectrum in SPH domain

%degwgtS1 = (1./dnum).^1; 
%spectS1  = coswgtS.*degwgtS1;               % Test signal 1/r spectrum in SPH domain
%degwgtS2 = (1./dnum).^2; 
%spectS2  = coswgtS.*degwgtS2;               % Test signal 1/r^2 spectrum in SPH domain

ctheta  = vecs(3,:)'./xrad';
phi     = atan2(vecs(2,:), vecs(1,:));
maxpnd  = ceil(maxdegp/2) + 1;
maxpnd2 = ceil(maxdegp/4) + 1;

if ppars
    if optimize   % no input weights are given
        fprintf(1,'\n  Optimizing SPH fiters')
        fprintf(1,'\n  -----------------------------------------------------------')
        fprintf(1,'\n  Nvec  = %3i      maxdeg =%3i      degref =%5.1f      d0 =%2i' ,Nvec1,maxdeg,degref,d0)
        fprintf(1,'\n  Nlost = %3i      cosexp =%5.2f       wrW =%5.2f \n\n'         ,Nlost,cosexp,wrW)
    else
        fprintf(1,'\n  Evaluating SPH fiters')
        fprintf(1,'\n  -----------------------------------------------------------')
        fprintf(1,'\n  Nvec  = %3i      maxdeg  =%3i'                                ,Nvec1,maxdeg)
        fprintf(1,'\n  Nlost = %3i      cosexpS =%5.2f       wrS =%5.2f \n\n'        ,Nlost,cosexpS,wrS)
    end
end

%%%% Compute spherical harmonics sample values for each input vector
B       = [];
W_      = [];
S_      = [];
dmark   = [];

for degree = degrees      % l in standard SPH notaion
    ndeg  = degree/2+1;
    pall  = legendre(degree,ctheta,'norm');
    order = [1:degree]';
    pord  = pall(order+1,:);
    SPHd0 = sqrt(2)*cat(3, cos(order.*phi).*pord, sin(order.*phi).*pord);
    SPHd0 = [pall(1,:); reshape(permute(SPHd0,[3 1 2]),2*degree,Nvec1)];
    B     = cat(2,B,SPHd0');                               % SPH basis vectors are columns
    W_    = cat(1,W_,ones(size(SPHd0,1),1)*sphwgt(ndeg));  % weights used for optimizing weights for each x
    S_    = cat(1,S_,ones(size(SPHd0,1),1)*spectS(ndeg));  % signal SPH spectrum amplitude
    dmark = cat(1,dmark,ones(size(SPHd0,1),1)*2*mod(degree/2,2) - 1);
end
C    = diag(cert);
B    = C*B;                   % multiply with measurement certainties
B    = B./sqrt(sum(B.^2,1));  % normalize SPH basis
BB   = B'*B;
BB0  = BB - eye(size(BB));    % zero diagonal correlaion matrix
W_   = W_ + 1e-6*max(W_);
W_   = W_/sqrt(mean(W_.^2));
W_2  = W_.^2;
W    = diag(W_);
W2   = diag(W_2);
S_   = S_ + 1e-6*max(S_);
S_   = S_/sqrt(mean(S_.^2));
S    = diag(S_);
A    = B*W2*B';      % + diag((1 - cert)./(cert.^2 + 0.1));
A    = A + 1e-14*max(diag(A))*eye(size(A));  % add 'eps' to avoid singular matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Optimize filters, i.e. get dual SPH basis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if optimize  % no input weights are given
    B_dual   = A\(B*W2);     % weight vector set
else
    sndeg = 0;
    for degree = degrees
        ndeg   = degree/2+1;
        linpos = sndeg + [1:4*(ndeg-1)+1]';
        sndeg  = linpos(end);
        B_dual(:,linpos) = Vwgts0{ndeg};
    end
end
B_dual  = B_dual./sqrt(sum(B_dual.^2,1));  % normalize dual basis

%%%% store dual basis in cells
if optimize        % no input weights are given
    sndeg = 0;
    for degree = degrees
        ndeg   = degree/2+1;
        linpos = sndeg + [1:4*(ndeg-1)+1]';
        sndeg  = linpos(end);
        Vwgts{ndeg,1} = B_dual(:,linpos);
    end
else
    Vwgts = Vwgts0;
end

%%%% %%%% Get original and optimized SPH error for degree d
BtBd   = B'*B_dual;               % projection of SPH basis on weight vectors
BtBd   = BtBd./sqrt(diag(BtBd))./sqrt(diag(BtBd))'; % Normalize projections
BtBd0  = BtBd.*(1 - eye(size(BtBd)));               % zero diagonal correlaion matrix

wBB    = W*BB;
wBB0   = W*BB0;
sBB0   = S*BB0;

wBtBd  = W*BtBd;                                    % weighted signal responces
wBtBd0 = W*BtBd0;                                   % weighted error signal responces
sBtBd0 = S*BtBd0;

Bdcorr = B_dual'*B_dual;

for degree1 = degrees
    ndeg1   = degree1/2+1;
    linpos1 = [(degree1-1)*(degree1)/2 + 1 : (degree1+1)*(degree1+2)/2]'; % all SPH of degree1
    %lpos1   = linpos1'
    for degree2 = degrees
        ndeg2   = degree2/2+1;
        linpos2 = [(degree2-1)*(degree2)/2 + 1 : (degree2+1)*(degree2+2)/2]'; % all SPHs of degree1
        %lpos2   = linpos2'
        
        % optimized SPH errors
        B1tBd2   = sBtBd0(linpos1,linpos2);
        E_B1tBd2 = sum(sum(abs(B1tBd2).^2));       % sum errors for all orders for given degrees
        Err_SPH_1(ndeg1,ndeg2) = sqrt(E_B1tBd2);
        ErrSPH{ndeg2,1}(ndeg1) = Err_SPH_1(ndeg1,ndeg2); % cell{:,1} holds optimized errors
        
        % Original errors
        B1B2    = sBB0(linpos1,linpos2);
        E_B1B2  = sum(sum(abs(B1B2).^2));
        Err_orig_1(ndeg1,ndeg2) = sqrt(E_B1B2);
        ErrSPH{ndeg2,2}(ndeg1)  = Err_orig_1(ndeg1,ndeg2);  % cell{:,2} holds orig errors
        
        % Weighted optimized SPH errors
        wB1tBd2   = wBtBd0(linpos1,linpos2);
        wE_B1tBd2 = sum(sum(abs(wB1tBd2).^2));       % sum errors for all orders for given degrees
        wErr_SPH_1(ndeg1,ndeg2) = sqrt(wE_B1tBd2);
        ErrSPH{ndeg2,3}(ndeg1)  = wErr_SPH_1(ndeg1,ndeg2);  % cell{:,2} holds orig errors

        % Weighted original errors
        wB1B2    = wBB0(linpos1,linpos2);
        wE_B1B2  = sum(sum(abs(wB1B2).^2));
        wErr_orig_1(ndeg1,ndeg2) = sqrt(wE_B1B2);
        ErrSPH{ndeg2,4}(ndeg1)   = wErr_orig_1(ndeg1,ndeg2);  % cell{:,2} holds orig errors
    end
end

% Get cumulative errors
Err_SPH   = sqrt(cumsum(Err_SPH_1.^2,1));
Err_orig  = sqrt(cumsum(Err_orig_1.^2,1));
wErr_SPH  = sqrt(cumsum(wErr_SPH_1.^2,1));
wErr_orig = sqrt(cumsum(wErr_orig_1.^2,1));

Err_orig_1(1,1)  = nan;  % dc error for dc not meaningfull
wErr_orig_1(1,1) = nan;  % dc error for dc not meaningfull
Err_SPH_1(1,1)   = nan;  % dc error for dc not meaningfull
wErr_SPH_1(1,1)  = nan;  % dc error for dc not meaningfull
ErrSPH{1,1}(1,1) = nan;  % dc error for dc not meaningfull
ErrSPH{1,2}(1,1) = nan;  % dc error for dc not meaningfull
ErrSPH{1,3}(1,1) = nan;  % dc error for dc not meaningfull
ErrSPH{1,4}(1,1) = nan;  % dc error for dc not meaningfull
Err_orig(1,1)    = nan;  % dc error for dc not meaningfull
wErr_orig(1,1)   = nan;  % dc error for dc not meaningfull
Err_SPH(1,1)     = nan;  % dc error for dc not meaningfull
wErr_SPH(1,1)    = nan;  % dc error for dc not meaningfull

if peval > 0
    fprintf('\n  Total error for order:             0         2         4')
    fprintf('\n  -----------------------------------------------------------')
    fprintf('\n  Using weights from optimizer    %7.5f   %7.5f   %7.5f', wErr_SPH(end,1), wErr_SPH(end,2), wErr_SPH(end,3))
    fprintf('\n  Using "white" signal weights    %7.5f   %7.5f   %7.5f \n\n', Err_SPH(end,1), Err_SPH(end,2), Err_SPH(end,3))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Plot section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotlev > 0
    spmg = 0.12;
    spsz = 0.86;
    max995 = 1e2;
    min995 = 1e-8;
    maxErr = max([Err_orig(:); Err_SPH(:); wErr_orig(:); wErr_SPH(:)]); 
    minErr = min([Err_orig(:); Err_SPH(:); wErr_orig(:); wErr_SPH(:)]); 

    figure(fignr1)   % total errors for degree 0 to d
    subplotHK(1,2,1,spmg,spsz)
    semilogy(degrees(1:maxpnd),Err_SPH(1:maxpnd,1:maxpnd2),'linewidth',2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    semilogy(degrees(1:maxpnd),Err_orig(1:maxpnd,1:maxpnd2),'--','linewidth',1)
    semilogy(degrees(1:maxpnd),ones(maxpnd,1),'linewidth',0.6,'color',[0 0 0])
    xlim([-0.25 maxdegp+0.25])
    ylim([min995, max995])
    xticks(degrees(1:maxpnd))
    grid on
    Dtext = ones(size(degrees(1:maxpnd2)'))*'Filter deg';
    hleg995 = legend([Dtext num2str(degrees(1:maxpnd2)')],'location','southeast');
    set(hleg995,'fontsize',10)
    xlabel('Max signal SPH degree')
    title('Total "white" signal error vs signal max SPH degree')
    hold off

    subplotHK(1,2,2,spmg,spsz)
    semilogy(degrees(1:maxpnd),wErr_SPH(1:maxpnd,1:maxpnd2),'linewidth',2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    semilogy(degrees(1:maxpnd),wErr_orig(1:maxpnd,1:maxpnd2),'--','linewidth',1)
    semilogy(degrees(1:maxpnd),ones(maxpnd,1),'linewidth',0.6,'color',[0 0 0])
    xlim([-0.25 maxdegp+0.25])
    ylim([min995, max995])
    xticks(degrees(1:maxpnd))
    grid on
    Dtext = ones(size(degrees(1:maxpnd2)'))*'Filter deg';
    hleg995 = legend([Dtext num2str(degrees(1:maxpnd2)')],'location','southeast');
    set(hleg995,'fontsize',10)
    xlabel('Max Signal SPH degree')
    title('Total weighted error vs signal max SPH degree')
    hold off
end
if plotlev > 1
    figure(fignr1+1)    % errors for degree d
    subplotHK(1,2,1,spmg,spsz)
    semilogy(degrees(1:maxpnd),Err_SPH_1(1:maxpnd,1:maxpnd2),'linewidth',2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    semilogy(degrees(1:maxpnd),Err_orig_1(1:maxpnd,1:maxpnd2),'--','linewidth',1)
    semilogy(degrees(1:maxpnd),ones(maxpnd,1),'linewidth',0.6,'color',[0 0 0])
    xlim([-0.25 maxdegp+0.25])
    ylim([min995, max995])
    xticks(degrees(1:maxpnd))
    grid on
    Dtext = ones(size(degrees(1:maxpnd2)'))*'Filter deg';
    hleg995 = legend([Dtext num2str(degrees(1:maxpnd2)')],'location','southeast');
    set(hleg995,'fontsize',10)
    xlabel('Signal SPH degree')
    title('Error vs signal SPH degree')
    hold off

    subplotHK(1,2,2,spmg,spsz)
    semilogy(degrees(1:maxpnd),wErr_SPH_1(1:maxpnd,1:maxpnd2),'linewidth',2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    semilogy(degrees(1:maxpnd),wErr_orig_1(1:maxpnd,1:maxpnd2),'--','linewidth',1)
    semilogy(degrees(1:maxpnd),ones(maxpnd,1),'linewidth',0.6,'color',[0 0 0])
    xlim([-0.25 maxdegp+0.25])
    ylim([min995, max995])
    xticks(degrees(1:maxpnd))
    grid on
    Dtext = ones(size(degrees(1:maxpnd2)'))*'Filter deg';
    hleg995 = legend([Dtext num2str(degrees(1:maxpnd2)')],'location','southeast');
    set(hleg995,'fontsize',10)
    xlabel('Signal SPH degree')
    title('Weighted error vs signal SPH degree')
    hold off

    figure(fignr1+2)
    %plot([1:size(W_)],[W_ ones(size(W))],'linewidth',1.5)
    semilogy([1:size(W_)],[W_ S_ ones(size(W_)) ones(size(W_))],'linewidth',1.5)
    xlim([0.5 size(W_,1)+0.5])
    ylim([1e-4 10])
    grid on
    xlabel ('SPH basis number')
    title('Weighting function (Expected signal strength)','fontsize',13)

    figure(fignr1+3)
    npldeg  = 2*(maxpnd-1);
    sphnum  = (npldeg-1)*npldeg/2;
    plBB    = BB(1:sphnum,1:sphnum);
    plwBB   = wBB(1:sphnum,1:sphnum);
    plBtBd  = BtBd(1:sphnum,1:sphnum);
    plwBtBd = wBtBd(1:sphnum,1:sphnum);
    max988  = max(abs([plBB(:); plBtBd(:)]));
    max988w = max(abs([plwBB(:); plwBtBd(:)]));

    spsz = 0.9;
    spd  = 0.3;
    subplotHK(4,1,1,spd,spsz)
    plot(plBB)
    xlim([0.5 size(plBB,1)+0.5])
    ylim(1.08*[-max988 max988])
    title('BB')

    subplotHK(4,1,2,spd,spsz)
    plot(plBtBd)
    xlim([0.5 size(plBB,1)+0.5])
    ylim(1.08*[-max988 max988])
    title('BtBd')

    subplotHK(4,1,3,spd,spsz)
    plot(plwBB)
    xlim([0.5 size(plBB,1)+0.5])
    ylim(1.08*[-max988w max988w])
    title('wBB')

    subplotHK(4,1,4,spd,spsz)
    plot(plwBtBd)
    xlim([0.5 size(plBB,1)+0.5])
    ylim(1.08*[-max988w max988w])
    title('wBtBd')

    figure(fignr1+4)
    subplotHK(1,2,1)
    imagesc(Err_orig.^0.33)
    colormap gray
    axis image

    subplotHK(1,2,2)
    imagesc(Err_SPH.^0.33)
    colormap gray
    axis image

    figure(fignr1+5)
    DM   = dmark*dmark';
    ccol = 1.3;
    subplotHK(1,2,1)
    gopimage(exp(i*ccol*DM),1.15)
    axis image

    subplotHK(1,2,2)
    RGB = Gop_RGB_func(abs(BtBd0).*exp(i*ccol*DM),[],0.2,0.25);
    image(RGB)
    axis image

    figure(fignr1+6)
    subplotHK(1,2,1)
    gopimage(Bdcorr.^2./abs(Bdcorr).^1.5)
    axis image

    subplotHK(1,2,2)
    gopimage(BtBd.^2./abs(BtBd).^1.5)
    axis image
    h989=figure(fignr1+7);
    %get(h989)
    pplBB    = interp2(abs(plBB),2,'nearest');
    pplwBB   = interp2(abs(plwBB),2,'nearest');
    pplBtBd  = interp2(abs(plBtBd),2,'nearest');
    pplwBtBd = interp2(abs(plwBtBd),2,'nearest');
    ppx      = interpn([1:size(plBB,2)],2,'linear');
    ppy      = interpn([1:size(plBB,1)],2,'linear');
    lip1     = [0 -10 1];
    lip2     = [-70 -70 1];
    lip3     = [50 -200 1];
    lip4     = [-100 -100 1];
    vang     = [-30,50];
    ast      = 0.1;
    dst      = min(1,0.3 + 0.0025*max(ppx));  % (don't get why this is needed)
    sst      = 0.0;
    colexp   = 0.12;
    magexp   = 0.33;
    max989   = 1.01*max988^magexp;
    xylim989 = [0.5 size(plBB,1)+0.5];
    calf     = -3.4;
    coff     = 2.2;
    csat     = 0.8;
    gim0     = 1e-5;

    subplotHK(2,2,1)
    gim = abs(pplBB);
    RGB = Gop_RGB_func(exp(i*(coff+calf*(gim+gim0).^colexp)),[],csat);
    surf(ppx,ppy,gim.^magexp,RGB,'ambientstrength',ast,'diffusestrength',dst,'specularstrength',sst,'edgealpha',0)
    shading('interp')
    xlim(xylim989)
    ylim(xylim989)
    zlim([0 max989])
    view(vang)
    light('position',lip1)
    light('position',lip2)
    light('position',lip3)
    %light('position',lip4)
    ylabel('Signal SPH number')
    xlabel('SPH filter number')
    title(['ProjectionMag^{',num2str(magexp),'}  (White signal, Origial coeffs)'],'fontsize',13)

    subplotHK(2,2,2)
    gim = abs(pplBtBd);
    RGB = Gop_RGB_func(exp(i*(coff+calf*(gim+gim0).^colexp)),[],csat);
    surf(ppx,ppy,gim.^magexp,RGB,'ambientstrength',ast,'diffusestrength',dst,'specularstrength',sst)
    shading('interp')
    xlim(xylim989)
    ylim(xylim989)
    zlim([0 max989])
    view(vang)
    ylabel('Signal SPH number')
    xlabel('SPH filter number')
    light('position',lip1)
    light('position',lip2)
    light('position',lip3)
    %light('position',lip4)
    title(['ProjectionMag^{',num2str(magexp),'}  (White signal, Optimized coeffs)'],'fontsize',13)

    h=subplotHK(2,2,3);
    gim = abs(pplwBB);
    RGB = Gop_RGB_func(exp(i*(coff+calf*(gim+gim0).^colexp)),[],csat);
    surf(ppx,ppy,gim.^magexp,RGB,'ambientstrength',ast,'diffusestrength',dst,'specularstrength',sst)
    shading('interp')
    xlim(xylim989)
    ylim(xylim989)
    zlim([0 max989])
    view(vang)
    ylabel('Signal SPH number')
    xlabel('SPH filter number')
    set(h,'Clipping','off')
    light('position',lip1)
    light('position',lip2)
    light('position',lip3)
    %light('position',lip4)
    title(['ProjectionMag^{',num2str(magexp),'}  (Weighted signal, Origial coeffs)'],'fontsize',13)

    h=subplotHK(2,2,4);
    gim = abs(pplwBtBd);
    RGB = Gop_RGB_func(exp(i*(coff+calf*(gim+gim0).^colexp)),[],csat);
    surf(ppx,ppy,gim.^magexp,RGB,'ambientstrength',ast,'diffusestrength',dst,'specularstrength',sst)
    shading('interp')
    xlim(xylim989)
    ylim(xylim989)
    zlim([0 max989])
    view(vang)
    ylabel('Signal SPH number')
    xlabel('SPH filter number')
    set(h,'Clipping','off')
    light('position',lip1)
    light('position',lip2)
    light('position',lip3)
    %light('position',lip4)
    title(['ProjectionMag^{',num2str(magexp),'}  (Weighted signal, Optimized coeffs)'],'fontsize',13)

end

if plotlev > 2
    figure(fignr1+91); mygopimage(W,0.7,0.25,2)
    figure(fignr1+92); mygopimage(BB,0.7,0.25,2)
    figure(fignr1+93); mygopimage(BtBd,0.7,0.25,2)
    figure(fignr1+94); mygopimage(wBB,0.7,0.25,2)
    figure(fignr1+95); mygopimage(wBtBd,0.7,0.25,2)
    figure(fignr1+97); mygopimage(wBB0,0.7,0.25,2)
    figure(fignr1+98); mygopimage(wBtBd0,0.7,0.25,2)
    figure(fignr1+99); mygopimage(sBtBd0,0.7,0.25,2)
end

if rendSPH > 0 | plotlev > 2
    alfa1   = 0.9;
    alfa2   = 0.05;
    bsz     = 0.017;
    rballsz = min(12*(240/Nvec1)^(1/3),29);
    shexp2  = 2;
    sballs  = 1;
    lincol  = 0;
    metric  = [8 1 1 1 2];
    chdens  = [1.00 1.50 0.03 0];
    cplot   = 2*cert - 1;                    % to get white for missing data
    relcol  = 1;
    Npdeg   = max(rendSPH,1);
    
    figure(fignr1+8)
    nB      = 0;
    x_2     = [vecs -vecs];
    cplot_2 = [cplot cplot];
    for ndeg = 1:Npdeg
        SPHdeg = [Vwgts{ndeg}; Vwgts{ndeg}];
        SPHnrm = SPHdeg/(max(abs(SPHdeg(:)))+eps);
        SPHnrm = (SPHnrm + 1)/2;           
        %nsp    = 4*(Npdeg-1)*(ndeg-1) + 2*(Npdeg-1)+1 - ndeg;  % symmetric pyramid
        nsp    = (4*(Npdeg-1)+1)*(ndeg-1);                    % left to right for each degree
        for ord = 1:4*(ndeg-1)+1;
            nB  = nB + 1;
            nsp = nsp+1;
            subplotHK(Npdeg,4*(Npdeg-1)+1,nsp)
            ballcolopt = cplot_2.*SPHnrm(:,ord)';  % to check B (should be the same as SPHdeg(ord,:)
            
            plot_3D_distribution(real(x_2)','colval',ballcolopt,'ashell',[alfa1 alfa2],'rshell',[0.99 1-bsz],'rhalf',1e9,'ballsize',rballsz,...
                                      'shellexp',[2 shexp2],'relcolscale',relcol,'shinyballs',sballs,'linecol',lincol,'metric',metric,'chdens',chdens);
            %plotellipticsamples_5plus3(real(x_2)',ballcolopt,[alfa1 alfa2],[0.99 1-bsz],...
            %                           1e9,rballsz,[2 shexp2],0,sballs,lincol,0,metric,chdens);
            %colorbar
            drawnow
        end
    end
end

if saveCode
    %save the code generating the results
    runnum = 88;
    fName  = mfilename;
    rName  = ['_',num2str(runnum)];
    fNameD = [fName,rName,'.data'];
    save(fNameD,'vecs','Vwgts');
    command = ['cp ' [fName,'.m'] ' ' [fName,rName,'.m']];
    system(command);
end


return














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                          % all terms for each unique monomial term.
    clear oprt
    clear oprt0
    clear oprtw
    Rexp = rexp*ones(size(xexp));
    for m = 1:size(x,2) % nb - exponents xexp,... are column vectors
        oprt0(:,m) = wt.*xrad(1,m).^rexp.*xhat(1,m).^xexp.*xhat(2,m).^yexp.*xhat(3,m).^zexp;
        if numel(Vwgts) > 1
            oprt(:,m) = bsxfun(@times, Vwgts(m), oprt0(:,m));  % Multiply with weights
        else
            oprt(:,m) = oprt0(:,m);
        end
        oprtw(:,m) = xrad(1,m).^Rexp;  % radial weights
        %oprt(:,m) = wt.*xhat(1,m).^xexp.*xhat(2,m).^yexp.*xhat(3,m).^zexp;
    end
    B = cat(1,B,oprt);% SPH transform matrix for the given x orientations
    W = cat(1,W,1./(1+(order/2)^wrW)*ones(nterms,1)); % weights to be used for optimizing weights for each x
    oprc  = sum(oprt,2)./(mean(oprtw,2));             % radially weighted linear sum of all identical terms
    oprE0 = sum(oprc(:).^2)/size(x,2);                % total energy for each outer order
    oprE(1+order/2) = oprE0*(order+1)/2;              % corrected to have the same 'DC' for each order

E_n = [oprE(1) diff(oprE)];                                     % Energy in each spherical harmonic shell
E_n = abs(E_n)/(1+1e-10) + 1e-12.*[0:maxorder/2]/2;   % add approximate uncertainty of computation (+ abs!!!)
S_n = sqrt(E_n);                                      % signal strength

owgt1    = coswgt.*(2 + norder.^2).^0;;
w_err(1) = sqrt(sum(E_n.*owgt1)/sum(owgt1));

owgt2    = coswgt./(2 + norder.^2).^1;
w_err(2) = sqrt(sum(E_n.*owgt2)/sum(owgt2));

owgt3    = coswgt./(2 + norder.^2).^2;
w_err(3) = sqrt(sum(E_n.*owgt3)/sum(owgt3));

owgt4    = coswgt./(2 + norder.^2).^3;
w_err(4) = sqrt(sum(E_n.*owgt4)/sum(owgt4));

Rotinv = real(exp(-14*S_n));  % makes 0.05 -> 1/2, resonable?
%rotinv_range = 2*sum(Rotinv);
%fprintf('\n Approximate  range  of  rotational invariance = %1.4f ', rotinv_range)

%%% prints and plots

if fignr < 6
    holdfig = 1;
    fignr   = abs(fignr);
else
    holdfig = 0;
end

if peval > 0
    fprintf('\n Average rotational sensitivity (1/n^0 weight)  = %1.2e', w_err(1))
    fprintf('\n Average rotational sensitivity (1/n^1 weight)  = %1.2e', w_err(2))
    fprintf('\n Average rotational sensitivity (1/n^2 weight)  = %1.2e', w_err(3))
    fprintf('\n Average rotational sensitivity (1/n^3 weight)  = %1.2e', w_err(4))
    
    fprintf('\n Range of rotational invariance  30%% tolerance  =%5.1f ', findtol(S_n,0.3))
    fprintf('\n Range of rotational invariance  10%% tolerance  =%5.1f ', findtol(S_n,1e-1))
    fprintf('\n Range of rotational invariance   1%% tolerance  =%5.1f ', findtol(S_n,1e-2))
    fprintf('\n Range of rotational invariance  e-3 tolerance  =%5.1f ',  findtol(S_n,1e-3))
    fprintf('\n Range of rotational invariance  e-4 tolerance  =%5.1f ',  findtol(S_n,1e-4))

    if maxorder > 3
        fprintf('\n 4:th order rotational sensitivity (ref cnd #)  = %1.2e',      S_n(2))
    end
    fprintf('\n 2:th order rotational sensitivity (ref rotinv) = %1.2e \n\n', S_n(1))

    if 0 == 1
        figure(fignr)  % plot rotinv error amplitude for spherical harmonic orders 0-maxorder
        subplot(2,1,1)
        plot([2:2:maxorder],log10(max(eps,[S_n])))%,'color',color)
        xlim([0 maxorder])
        ylim([-6 0])
        grid on
        title(' log10 of relative rotational rms error for each order ','fontsize',16,'fontweight','normal')
        if holdfig
            hold on
        else
            hold off
        end

        subplot(2,1,2)
        % plot rotational invariance estimate agains order
        plot([0:2:maxorder],[1 Rotinv])%,'color',color)
        xlim([0 maxorder])
        ylim([0 1])
        title(' Special rotational invariance function','fontsize',16,'fontweight','normal')
        grid on
        %figtitle(' Rotational evaluation plots ')
        if holdfig
            hold on
        else
            hold off
        end
        drawnow
    end
end

%%%%%%%%%%%%%%%%%%% Optimize weights %%%%%%%%%%%%%%%%%%%%

if numel(Vwgts) == 1                      % only optimize when no weights are entered
    %q_ideal = [1; zeros(size(B,1)-1,1)];  % only 'DC' should have a conribution for perfect rotation invariance
    q_ideal = [1e-16; 0; 1; zeros(size(B,1)-3,1)];  % 2:nd order ideal
    %q_ideal = [1; 0.01*ones(5,1); zeros(size(B,1)-6,1)];  % 0 and 2:nd order
    W       = diag(W + 1e-6);
    Vwgts = pinv(B'*W*B)*(W*B)'*q_ideal;
    Vwgts  = Vwgts/sqrt(mean(Vwgts.^2));
end

function tolord = findtol(Sn,tol);
Snx = [0 Sn 1];
ix2 = min(find(Snx > tol));
ix1 = ix2 - 1;
dS1 = tol - Snx(ix1);
dS2 = Snx(ix2) - tol;
tolord = 2*((dS2*ix1 + dS1*ix2)/(dS1 + dS2 + eps) - 1);


return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ESTA #

R8XXR8R8R94FT72T
