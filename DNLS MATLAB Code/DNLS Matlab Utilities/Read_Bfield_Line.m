function [PStr, BN, MSG, flag] = Read_Bfield_Line(ArchiveBase, TrialNumber, RunNumber, LineNumber)
%[PStr1, B1, MSG1,flag1] = Read_Bfield_Line(ArchiveBase, TrialNumber1,RunNumber1, LineNumber1);
%   Given:  file address for Bfield (Bfile)
%           Bfield line number to read (LineN)
%           Parameters for trial & run (from PStr)
%           --> extract and return a complex valued array containing the B-field
%           --> create and return a message on the status of the operation
%
%   Error Checking:     check that the file exists
%                       check that the file contains the line number requested
%% To Do:   Work towards saving B-field files in the 'matfile' format
%           --> This allows just the necessary data to be loaded into
%           memory. The text file format requires the entirety of the file
%           to be loaded and processed. The more efficient memory management 
%           becomes more important given Matlab's (in conjunction with
%           Windows) failure to free memory once assigned to a cleared variable.
%           see: https://www.mathworks.com/help/matlab/ref/matfile.html
%% Construct File Addresses
TrialFolderRel = ['\Trial ',TrialNumber];%e.g. '\Trial 11'
TrialFolderAbs = [ArchiveBase,TrialFolderRel];%e.g. 'F:\ ...\Trial 11\'
Bfile = [TrialFolderAbs,'\Bfields\Bfield_',TrialNumber,'_',RunNumber,'.txt'];%e.g. 'F:\ ...\Trial 11\'  
BFieldFile = [TrialFolderAbs,'\Bfields\Bfield_',TrialNumber,'_',RunNumber,'.txt'];%e.g. 'F:\ ...\Trial 11\' 
BfieldMatfile = [TrialFolderAbs,'\Bfields\Bfield_',TrialNumber,'.mat'];%e.g. 'F:\ ...\Trial 11\'
thisFile = mfilename();

%***********************************************************    
 %Read Excel File with Trial and Run Parameters
%***********************************************************     
['From: ',thisFile,' Begin reading parameters ...  ']
    %Read Excel File with Trial and Run Parameters
        % Extract PStr: the structure with members such as Physics, Numerics etc.
        % with each member containing {keys, values} of parameters describing a given trial & run  
        if exist(BfieldMatfile,'file')
            ['From: ',thisFile,' Begin reading parameters from Matfile Archive  ']
            BFMat = matfile(BfieldMatfile);
            PStr = BFMat.PStr(str2num(RunNumber),1);
        else 
            ['From: ',thisFile,' Begin reading parameters from Excel Archive  ']
            Efile = [TrialFolderAbs,'\*.xlsx'];%wild card for files in Trial folder with '.xls' extension
            RelXLfile = dir(Efile);% returns, presumably, the single file with '.xlsx' extension
            XLfile = [TrialFolderAbs,'\',RelXLfile.name];%absolute address of excel file with trial parameters

            %Read sheet 'RunNumber' of the Excel file 'XLfile' and create
            %'PStr' ... a  Parameter structure defined as:
            % Pstr = struct('Dcon',Diagnostics, 'Fcon', Folder,'GCon',
            % Graphics, 'NCon',Numerics,'PCon',Physics,'TCon',Trial)
            % To find the value of 'L' stored in the Numerics container use:
            % L = Pstr.NCon('L') ... for example

            [PStr] = DNLS_Excel_Read_v3(XLfile, RunNumber);
        end


    if exist(BfieldMatfile,'file') && 0
        ['From: ',thisFile,' Begin Read B-field Line. Trying Matfile Archive  ']
        BFMat = matfile(BfieldMatfile);
        BfieldNames = fieldnames(BFMat);
        idxRun = strfind(BfieldNames,'Run');
        if str2num(RunNumber) <= length(BFMat.PStr(:,1))
            runN = ['Run',RunNumber];
        end
        if sum(~cellfun('isempty',strfind(BfieldNames,runN))) > 0
            if LineNumber <= length(BFMat.(runN)(:,1))              
                B4 = BFMat.(runN)(LineNumber,:);
            else
                ['Line # not in Matfile']
                length(BFMat.(runN)(:,1))
                return;
            end
        else
            ['Run # not in file/nr']
            ['Field Names:']
            fieldnames(BFMat)
            return;
        end
    else 
        %Read all lines from Bfield file for given TrialNumber and RunNumber
        ['From: ',thisFile,' Begin Read B-field Line. Trying text file archive']
            BFh = fopen(BFieldFile);
            if BFh < 1
                [BFieldFile, ' does not exist.']
                return
            end
            BFline = fgetl(BFh);
            BFlines = cell(0,1);
            NumOfLines = 0;
            while ischar(BFline)
%                 BFlines{end+1,1} = BFline;
                BFline = fgetl(BFh);
                if strfind(BFline,'LineStart')>0
                    NumOfLines = NumOfLines + 1;
                    BFlines{end+1,1} = fgetl(BFh);
%                     clc;
                    ['From: ',thisFile,' verifying B-field Line # ',num2str(NumOfLines)]
                end
            end
            fclose(BFh);

        %Find Bfield profile for Line n  
            LineStartN = ['LineStart',num2str(LineNumber)];
            numlines = length(BFlines); %returns number of lines in BFh
            IndexBF = strfind(BFlines, LineStartN); %returns line # of flag before Bfield Line n
            lineN = find(not(cellfun('isempty',IndexBF))) + 1; %gives line # of nth Bfield Line
            if LineNumber < (NumOfLines + 1)
%                 clc;
                ['From: ',thisFile,' loading B-field Line # ',num2str(LineNumber)]
                temp = BFlines(LineNumber);
                C = textscan(temp{1},'%f');
                B4 = C{1}.'; %B is an N/4 complex valued, double precision array
            else
                ['Line # ', num2str(LineNumber),' is greater than the ', num2str(NumOfLines), ' lines in the file.' ]
                LineStartN = ['LineStart',num2str(numlines)];
                LineNumber = numlines;
                ['From: ',thisFile,' loading B-field Line # ',num2str(LineNumber)]
                temp = BFlines(LineNumber);
                C = textscan(temp{1},'%f');
                B4 = C{1}.'; %B is an N/4 complex valued, double precision array
            end
    end

%% Try to Open Bfile
    BfieldID = fopen(Bfile);
    BFlines = cell(0,1); % Declare BFlines as null cell array
    BN = -1;
    flag = -3;
   
    if BfieldID > 0
         % File exists.  Do stuff....
      MSG = sprintf('File exists:\n%s', Bfile);
      flag = 1;
    else
      % File does not exist.
      MSG = sprintf('Warning: file does not exist:\n%s', Bfile);
      flag = -1;
      return;
    end
%% Load all lines of Bfile into BFlines cell array & close Bfile    
    BFline = fgetl(BfieldID);
%     BFlines = cell(0,1);
    while ischar(BFline)
        BFlines{end+1,1} = BFline;
        BFline = fgetl(BfieldID);
    end
    fclose(BfieldID);
%% Check if LineN is found in BFlines        
    LineStartN = ['LineStart',num2str(LineNumber)];
    numlines = length(BFlines); %returns number of lines in BFh
    IndexBF = strfind(BFlines, LineStartN); %returns line # of flag before Bfield Line n
    lineN = find(not(cellfun('isempty',IndexBF))) + 1; %gives line # of nth Bfield Line
    
    if isempty(lineN) == 1
        MSG = [MSG, sprintf('\nBfield line # %s NOT found.\n', num2str(LineNumber))];
        flag = -2;
        return
    else
        MSG = [MSG, sprintf('\nBfield line # %s found.\n', num2str(LineNumber))];
    end
%% Save BN = B(lineN) as a complex valued array of doubles    
    Btemp = BFlines(lineN);
    CN = textscan(Btemp{1},'%f');
    BN = CN{1}.'; %B is an N/4 complex valued, double precision array
end