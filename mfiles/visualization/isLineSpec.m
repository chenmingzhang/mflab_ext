function [LS,c,m,L] = isLineSpec(lineSpec)
%ISLINESPEC check if lineSpec is a legal Matlab lineSpec
%
% Example:
%    [LS,c,l,m] = isLineSpec(lineSpec) -- 
%
%    legal lineSpecs are a combination of
%    a color from       'bwgkmcyw'
%    a line type from   '--',':','-.','-','none'
%    a marker from      '.os^vp+x*'
%
% see also: mf_color mf_lineType isColor
%
% TO 210101

if ~ischar(lineSpec)
    LS = false;
    c='';
    m='';
    L={};
    return;
end

LSPECS  = {'--',':','-.','-'};
COLORS  = 'brgkmcyw';
MARKERS = '.os^vp+x*';

j = strfind(lineSpec,'none');
if ~isempty(j)
    lineSpec(j:j+3)=[];
end

%% Does lineSpec have colors ?
clr      = ismember(lineSpec,COLORS);
c        = lineSpec( clr);
lineSpec = lineSpec(~clr);

%% Does lineSpec have markers ?
mrk      = ismember(lineSpec,MARKERS);
m        = lineSpec( mrk);
lineSpec = lineSpec(~mrk);

%% Does lineSpec have lineStyles ?
L{numel(LSPECS),1} = [];
for i=1:numel(LSPECS)
    j=strfind(lineSpec,LSPECS{i});
    if ~isempty(j)
        L{i} = LSPECS{i};
        lineSpec(j:j+length(LSPECS{i})-1)=[];
    end
end
L = L(~cellfun(@isempty,L));

LS = isempty(lineSpec);

if nargout>1
    if ~isempty(c), c=c(1); end
    if ~isempty(m), m=m(1); end
    if ~isempty(L), L= L{1}; end
end