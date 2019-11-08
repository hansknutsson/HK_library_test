function [Heqcols] = Hist_eq_cols(cols,varargin)
% --------
% function [Heqcols] = Hist_eq_cols(cols,varargin)
% --------
% Makes 'soft' histogram equalization to spread values more evenly
%
% @ H Knutsson   Jan 2019
%
% Obtional arguments:XD
% Histsz  = 200;
% Hfiltsz = 15;
% Histexp = 0.25;

Histsz  = 200;
Hfiltsz = 15;
Histexp = 0.25;
for nvarin = 1:2:numel(varargin)
    switch varargin{nvarin}
      case 'Histsz';       Histsz   = varargin{nvin+1};
      case 'Hfiltsz';      Hfiltsz  = varargin{nvin+1};
      case 'Histexp';      Histexp  = varargin{nvin+1};
        error(['Unexpected option: ' varargin{nvin}])
    end
end
colHist = (conv(hist(cols,Histsz),ones(1,Hfiltsz),'full')./...
           conv(ones(size(hist(cols,Histsz))),ones(1,Hfiltsz),'full')).^Histexp;
cumH    = (cumsum(colHist) - colHist(1))/(sum(colHist) - colHist(1));
coltab  = (cumH + eps)/(1+eps);
Heqcols = coltab(1+round((size(coltab,2)-1)*(cols-min(cols)+eps)/...
                           (max(cols)-min(cols)+eps)))+eps;
