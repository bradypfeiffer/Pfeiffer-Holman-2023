function [PStr   , BN        , LineNselected    , MSG , flag ] = Continuation_Read_Bfield_Line(ArchiveBaseOld, ArchiveBaseNew, TrialNumberOld, RunNumberOld, LineNumber   , TrialNumberNew, RunNumberNew, Verbal    )
%%        [PStrOld, B4NewStart, LineNumberOldRead, MSG1, flag1] = Continuation_Read_Bfield_Line(ArchiveBase, TrialNumberOld, RunNumberOld, LineNumberOld, TrialNumber   , RunNumber_start );
%   Given:  ArchiveBase = absolute file address for '...\Archive\Bobs_Trials' containing all archived trial runs to be referenced
%
%           TrialNumberOld = '16' ... the trial number as a char string for archived trial to read
%           RunNumberOld   = '10' ... the run number as a char string for archived trial to read
%           LineNumber     =  34   ... the line number as a double to be read from the archived trial
%
%           TrialNumberNew = '16' ... the trial number as a char string for new trial to write
%           RunNumberNew   = '10' ... the run number as a char string for new run to write
%***************************************************************************
%     If TrialNubmerNew = 0 ... zero as a double
%     --> only LineNumberOld from TrialNumberOld and RunNumberOld 
%         will be read and returned and no lines will be written to a new Trial and Run
%***************************************************************************


%file address for Bfield (Bfile)
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
nl = newline;
TrialFolderRelOld = ['\Trial '  , TrialNumberOld   ]; %e.g. '\Trial 11'
TrialFolderAbsOld = [ArchiveBaseOld, TrialFolderRelOld]; %e.g. 'F:\ ...\Trial 11\'
BfileOld  = [TrialFolderAbsOld, '\Bfields\Bfield_', TrialNumberOld, '_', RunNumberOld, '.txt' ]; %e.g. 'F:\ ...\Trial 11\'  
BfieldMatfileOld = [TrialFolderAbsOld, '\Bfields\Bfield_', TrialNumberOld, '.mat' ]; %e.g. 'F:\ ...\Trial 11\'

% if TrialNumberNew == 0
    TrialFolderRelNew = ['\Trial '  , TrialNumberNew   ]; %e.g. '\Trial 11'
    TrialFolderAbsNew = [ArchiveBaseNew, TrialFolderRelNew]; %e.g. 'F:\ ...\Trial 11\'
    BFieldNew = [TrialFolderAbsNew, '\Bfields\Bfield_', TrialNumberNew, '_', RunNumberNew, '.txt' ]; %e.g. 'F:\ ...\Trial 11\'
    BfieldMatfileNew = [TrialFolderAbsNew, '\Bfields\Bfield_', TrialNumberNew, '.mat' ]; %e.g. 'F:\ ...\Trial 11\'
% end

thisFile = mfilename();

%***********************************************************    
 %Read Excel File with Trial and Run Parameters
%***********************************************************     
if Verbal
    ['From: ',thisFile,' Begin reading parameters ...  ']
end
    %Read Excel File with Trial and Run Parameters
        % Extract PStr: the structure with members such as Physics, Numerics etc.
        % with each member containing {keys, values} of parameters describing a given trial & run  
        BFMatWorks = 0; % Assuming BfieldMatfile does not work for some reason
        ExcelWorks = 0;
        runNOld = ['Run',RunNumberOld];
        MSG = "Whaaat?";
        
        % Check if Matfile exists
        if exist(BfieldMatfileOld,'file')            
            BFMat = matfile(BfieldMatfileOld);
            BFMatWorks = 1;
            varList = who(BFMat);
            if Verbal
                ['From: ',thisFile,' Begin reading parameters from Matfile Archive  ']
            end
            MSG = ['Archived Matfile exists'];
        else
            if Verbal
                ['File does not exist:',newline, BfieldMatfileOld]
            end
            BFMatWorks = 0;
            MSG = ['Archived Matfile does not exist.'];
        end
        
        % Check if RunN is in Matfile
        runNOld = ['Run',RunNumberOld];
        if BFMatWorks && ismember(runNOld,varList)
            BFMatWorks = 1;
            MSG = [MSG,', archived run # was found'];
            if Verbal
                [runNOld,' found in Matfile']
            end
        else
            if Verbal
                [runNOld,' is not found in',newline, BfieldMatfileOld]
            end
            BFMatWorks = 0;
            MSG = [MSG,' but archived run # was not found'];
        end
        
        % Check if PStr is in Matfile
        if ismember('PStr',varList)  && BFMatWorks
            PStr = BFMat.PStr(str2num(RunNumberOld),1);
            BFMatWorks = 1;
            if Verbal
                ['PStr found.'] 
            end
            MSG = [MSG,', PStr was found'];
        else
            if Verbal
                ['PStr not found in',newline,BfieldMatfileOld]
            end
            BFMatWorks = 0;
            MSG = [MSG,' PStr was not found'];
        end

    BN = -1;LineNselected = -1; flag = -1; 
    
    if (str2num(TrialNumberNew) ~= 0) && (BFMatWorks ~= 0) % If a new Trialnumber is given, write linenum = 1: (LineNumber -1) to new archive
        % Read in each line from archive 'RunNumber' and write to
        % to BMat and Bfile for 'RunNumberNew'
        for linenum = 1:(LineNumber - 1)
            Stemp1 = [nl,'LineStart',num2str(linenum),nl];
            BfileNew = fopen(BFieldNew,'a');
            fprintf(BfileNew,Stemp1);
            fclose(BfileNew);

            b4 = b(1:4:length(b)) + 1i*10^(-25);
            dlmwrite(BFfile,b4,'precision','%.6f','-append','delimiter',' ');
            BFMat.(['Run',RunNumberOld])(linenum,1:length(b4)) = b4;
        end  
        MSG = [MSG,' archived Bfield lines were written to Matfile for new trial'];
    else
        MSG = [MSG,' Bfield lines WERE NOT written to Matfile for new trial'];
    end
    
    if BFMatWorks % Evaluate and return BN: the archived LineNumber

        if LineNumber <= length(BFMat.(runNOld)(:,1))
            B4 = BFMat.(runNOld)(LineNumber,:);
            if Verbal
                ['Line # ',num2str(LineNumber),' was just read from BFMat.']
            end
            LineNselected = LineNumber;
            BN = B4;
        else
            LineNumber = length(BFMat.(runNOld)(:,1));
            if Verbal
                ['Line # ',num2str(LineNumber),' is the last line archived and this line will be used.']
            end
            LineNselected = LineNumber;
            B4 = BFMat.(runNOld)(LineNselected,:);
            BN = B4;
        end
        MSG = [MSG,' and Bfield line #',num2str(LineNselected),' was returned.'];
        flag = BFMatWorks;
    end

 
    
end