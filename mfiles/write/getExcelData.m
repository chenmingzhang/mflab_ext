function [Parnams,Parvals,Partxthdr,Partxt,singleValue]=getExcelData(XLSname,SHTname,direction,varargin)
%GETEXCELDATA reads info from excel worksheet
%
% Example:
%    [Parnams,Parvals,Partxthdr,Partxt,SingleValue]=getExcelData(XLSname,SHTname,direction[,fill|noFill][,variableName])
%
% What it does:
%   Get MODFLOW/MT3D/Seawat etc parameter names and parameter values stored
%   in the Microfost Excel worksheet given by workbook name XLSname, and
%   worksheet name SHTname.
%   direction is the direction in which the labels of the data are arranged.
%  'V[ertical]' means that the parametrs names are arranged vertically
%   with the values to the right of the parameter names so that multiple
%   columns are possible.
%   'H[orizontal]' means the other way around, that is, the parameter names
%   are arranged horizontally, with the variables themselves below the names,
%   allowing multiple rows, such as for layers and streess periods.
%
% varargin can be:
%   'fill' true    indicating that empty spreadsheet cells are filled using the last
%                  filled cell above. This is the default.
%   'fill' false   replaces empty cells with NaN
%
% variableName is optional name of single variable, whose value may be
% outputted as 5th output argument. Variable names are those in the
% corresponding worksheet.
%
% In vertical setup, lines with blank first cells, will be removed.
% This means that such lines are allowed on the sheet.
% Parnams are the names of the header
% Parvals are the numerical values
% ParText are the textual values
%
% TO 081231 110804 120415 130204 130626

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<3,
    error('getExelData:nargin:insufficientInput',...
        ['getExcelData: need at least three input arguments as follows from its usage:\n',...
        '[Parnams,Parvals,Partxthdr,Partxt,singleValue]=getExcelData(XLSname,SHTname,direction,variableName)\n']);
end

% Usefull if fill is on and NaN's are desired
NaNthreshold = 1e50;  % values for which abs(values)>NaNthreshold will be set to NaN

% option to fill or not fill open field in worksheet.
% default is fill.
[fill  ,  ~   ] =  getProp(varargin,'fill',true);

singleValueRequest = nargout>4 && nargin>3;

if singleValueRequest
    variableName = varargin{1};
    %varargin(4)=[];
end

%% Input checking to prevent trouble further down

% Third argument either H or V?
if isempty(regexp(direction,'^[hHvV]','once'))
   error('%s: The label direction (3rd argument) must be of type char and start with ''H'' or ''V''',mfilename);
end

%% first and second inputs must be of class char
if ~ischar(XLSname)
    error('%s: XLSname <<%s>> must be valid name of a Microsost Excel Spreadsheet',...
        mfilename,XLSname);
end
if ~ischar(SHTname)
    error('%s: SHTname <<%s>> must be a valid worksheetNm in workbook <<%s>>',...
        mfilename,XLSname,SHTname);
end

%% Second arg must be a worksheet in workbook given in first argument
[STATUS,SHEETS] = xlsfinfo(XLSname);
if isempty(regexp(STATUS,'Microsoft','once'))
    error('%s: XLSname %s not a legal %s',mfilename,'Microsoft Excel Spreadsheet');
end
if ~ismember(SHTname,SHEETS)
    error(['%s: SHTname %s not in workbook %s\n'...
           'REMEDY: maybe the workbook was not saved as an\n',...
           'Excel 5.0/95 workbook\n',...
           'The latter is necessary on non-Windows PC''s because\n',...
           'the Matlabs xlsread can''t read more recent fileformats of\n',...
           'Excel unless of version 5.0/95.\n',...
           'This will change depending on the Mathworks.'],...
           mfilename,SHTname,XLSname);
end

%% Ok, continue
direction = upper(direction(1));

warning('off','all');

if ismac
    warning('off','MATLAB:xlsread:Mode')
    try
    [~,~,Raw]=xlsread(XLSname,SHTname,'','basic');
        
    catch ME
   %     [STATUS,SHEETS] = xlsfinfo(XLSname);
        error(['%s: %s\n\n',...
            'REMEDY: First make sure the spelling of workbook and sheetName is exact.\n',...
            '   Notice that this spelling is case sensitive.\n',...
            '   Next, notice that this may happen on the MAC or other non-PC operating systems',...
            '   lacking excel.com so that xlsread must be run in basic mode.\n',...
            '   However, the Matlab function xlsread when running in basic mode\n',...
            '   will only work properly if the xls worksheet was saved\n',...
            '   as a <<Microsoft Excel 5.0/95 workbook>>.\n',...
            '   (On MS-Windows this (Matlab)-limitation does not apply).',...
            '   Therefore, saving the your Excel file in this format on Mac may solve this problem.\n',...
            '   You may use the command\n',...
            '   [~,sheets] = xlsfinfo(''%s'')\n',...
            '   to get a list of the sheets in the workbook that xlsread sees.\n',...
            '   My experience is that xlsread is not bug free, sometimes it does not see\n',...
            '   sheets that are definetly there. I have solved this by moving sheets around\n',...
            '   in the excel workbook (moving the invisible one to the front for instance,\n',...
            '   or copying the contents to a new sheet, remvoving the old one\n',...
            '   and renaming the new one to the old one.'],...
            mfilename,ME.message,XLSname);
    end
else
    try
        [~,~,Raw]=xlsread(XLSname,SHTname);
    catch ME
        error('%s: %s\nREMEDY: Check the exact spelling of the sheetName, this is case sensitive\n',...
            mfilename,ME.message);
    end
end

% Remove all empty columns
nans = cellfun(@(a) isnan(a(1)),Raw);
Raw  = Raw(~all(nans,2),:);

% Remove all empty columns
nans = cellfun(@(a) isnan(a(1)),Raw);
Raw  = Raw(:,~all(nans,1));


if direction=='H'

    txtCol  = all(cellfun( @(a) isnan(a(1)), Raw) | cellfun(@ischar,Raw),1);
    numCol  = ~txtCol;

    DataRow = all(cellfun(@isnumeric,Raw(:,numCol)),2);

    if ~any(DataRow)
        error(['%s: There are no numerical data rows in worksheet <<%s>> of workbook <<%s>>.\n',...
               'REMEDY: Make sure there is at least one numerical row in your worksheet.'],...               
            mfilename,SHTname,XLSname);
    end

    LblRow  = find(~DataRow);
    LblRow  = LblRow(end);

    if any(cellfun(@(a) isnan(a(1)), Raw(LblRow,:)))
        I = find(cellfun(@(a) isnan(a(1)),Raw(LblRow,:)));
        error(['There is at least one label that is EMPTY, or NUMERICAL or a FORMULA in\n', ...
               'column <<%d>> in your worksheet <<%s>> of workbook <<%s>>.\n',...
               'REMEDY: Verify that you don''t use character formulas in the labels in Excel (Matlab can''t read their textual values)'],...
                I(1),SHTname,XLSname);
    end

    Parnams    = Raw(LblRow(end),numCol);
    Parvals    = cell2mat(Raw(DataRow,numCol));
    Partxthdr  = Raw(LblRow(end),txtCol);
    Partxt     = Raw(DataRow    ,txtCol);
        
    %% Fill in all NaN cells as of row 2
    % TO 130612
    if fill
       if any(isnan(Parvals(1,:)))
            error(['When fill is on (default), Empty or NaN values are not allowed\n',...
                   'in the first row under the labels,\n',...
                   'in workbook <<%s>> sheet <<%s>> labels <<%s>>.\n',...
                   'Alternatively use value >%g to set cells to NaN''s'],...
                    XLSname,SHTname,sprintfs(' %s',Parnams(isnan(Parvals(1,:)))),NaNthreshold);
       end

        % fill empty numeric data lines if necessary
        for iL = 2:size(Parvals,1)
            I = isnan(Parvals(iL,:));
            if ~isempty(I)
                Parvals(iL,I)=Parvals(iL-1,I);
            end
        end

        % fill empty text lines if necessary
        if ~isempty(Partxthdr)
            for iL = 2:size(Partxt,1)        
                J = cellfun(@(a) isnan(a(1)), Partxt(iL,:));
                if ~isempty(J)
                    Partxt(iL,J) = Partxt(iL-1,J);
                end        
            end
        end
    end

    % potentially insert NaN by using an extreme threshold
    Parvals(abs(Parvals)>NaNthreshold)=NaN;

    if singleValueRequest
        I = strmatchi(variableName,Parnams);
        if I(1)
            singleValue = Parvals(:,I(1));
        elseif strmatchi(variableName,Partxthdr)
            singleValue = Partxt(:,strmatchi(variableName,Partxthdr));
        else
            error('Can''t find label <<%s>> in table in workbook <<%s>> worksheet <<%s>>.',...
                   variableName,XLSname,SHTname);
        end
    end
    
    % Cut off any trailing NaNs from Partxt
    lastNonNaN = find(~all(cellfun(@(a) isnan(a(1)) ,Partxt),2),1,'last');
    Partxt(lastNonNaN+1:end,:)=[];

else % direction == 'V'

    txtRow  = all(cellfun(@(a) isnan(a(1)),Raw) | cellfun(@ischar,Raw),2);
    numRow  = ~txtRow;
    if ~any(numRow)
            error('Make sure there is at least one numerical row in your data of workbook <<%s>> worksheet<<%s>>.',...
                   XLSname,SHTname);
    end

    % Determination of DataCol requires at least columns with at least one numeric value
    % one numeric row (after the label columns)
    DataCol = any(cellfun(@isnumeric,Raw(numRow,:)),1);
    LblCol  = ~DataCol; % columns with all labels (need the last)
    i=find(~LblCol,1,'first');
    if ~isempty(i) && i>1
        LblCol(i:end)=false;
    end

    if any(cellfun(@(a) isnan(a(1)),Raw(:,LblCol)))
        I = find(cellfun(@(a) isnan(a(1)) ,Raw(:,LblCol)));
        error(['There is at least one numerical label in row <<%d>> in your workbook <<%s>>, worksheet <<%s>>.\n',...
               'Verify that you don''t use character formulas in the labels, as Matlab can''t read their textual output'],...
                I(1),XLSname,SHTname);
    end

    Partxthdr  = Raw(txtRow,find(LblCol,1,'last'));
    Partxt     = Raw(txtRow,DataCol    );
    Parnams    = Raw(numRow,find(LblCol,1,'last'));

    % remove any cell with text from the numeric range
    for i=1:numel(numRow)
        for j=find(DataCol)
            if ischar(Raw{i,j})
                Raw{i,j}=NaN;
            end
        end
    end
    % then assign this now truly numeric range to Parvals
    Parvals    = cell2mat(Raw(numRow,DataCol));
    
    
    % Cut off any trailing NaNs from Partxt
    lastNonNaN = find(~all(cellfun(@(a) isnan(a(1)) ,Partxt),1),1,'last');
    Partxt(:,lastNonNaN+1:end)=[];

    if singleValueRequest
        I =strmatchi(variableName,Parnams);
        if I(1)
            singleValue = Parvals(I(1),:);
            singleValue(isnan(singleValue))=[];
        elseif strmatchi(variableName,Partxthdr)
            singleValue = Partxt(strmatchi(variableName,Partxthdr),:);
            singleValue(isnan(singleValue))=[];
        else
            error('Can''t find label <<%s>> in table in workbook <<%s>> worksheet <<%s>>.',...
                   variableName,XLSname,SHTname);
        end
    end
end

