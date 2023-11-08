%%Calls DNLS_Reader_vn(ArchiveBase, TrialNumber, RunNumber, LineNumber)
%           |-> DNLS_Excel_Read(XLfile, RunNumber)
%                   |-> xlsread(XLfile,Sheet,'D17:E22')
%********************************************************************
%% Changes in Call_DNLS_Reader v6 from v5:
%   
%   1) Create PStr within 'Call' Script
%       a) Build 'DNLS_Excel_Read_v3' call from DNLS_Reader_v5 into 'Call'
%       b) Move 'User Declared Parameters' from DNLS_Reader_v5 into 'Call'
%       c) Combine parameters of (a) and (b) into PStr parameter structure 
%           similar to format used in 'Call_DNLS_v3c'
%       d) Pass PStr from 'Call' to 'DNLS_Reader_vx'
%
%   2) Organize user flags (Video creation, Real EV, Complex EV, Newton's Method, Hybrid Box/Newton's Method)
%       for compatibility and formatting with PStr
%
%   3) Check for number of monitors to allow graphs to display on screen.
    %   Note that graphics information is derived from 'groot' through
    %   groot.MonitorPositions. Graphics settings are fixed upon Matlab startup
    %   and are fixed and, so, will not change if a monitor is added or
    %   removed. Matlab needs to be restarted to update groot. However, Matlab
    %   has published an API
    %   (https://www.mathworks.com/matlabcentral/fileexchange/31437-windowapi?s_tid=mwa_osa_a)
    %   that, evidently, can be used to check for system changes post-startup.
    %   It seems to require a C-compiler to be installed and accessible through
    %   Matlab for installation:
    %   (https://www.mathworks.com/help/matlab/matlab_external/choose-c-or-c-compilers.html)
%   
%% Changes in Call_DNLS_Reader v5 from v4:
%   Extends scattering data evaluation to the case of fully asymmetric
%   boundary conditions. Coding for this is in call to DNLS_Reader_v7
%% Changes in Call_DNLS_Reader v4 from v3:
%   1) Calls 'DNLS_Reader_v6':
%   Ability to read PStr parameters and B-field array from a '*.mat' mat
%   file format. This file format has the advantage of being a native
%   matlab file format and, as such, will function within either a Mac or
%   PC OS. Also, the mat file format allows a native, streamlined means of
%   assigning and accessing values from the file. Finally, it allows one
%   variable (or B-field line) to be read and stored at a time which should
%   reduce read/write time and free memory.
%% ToDo:
%
%   1)  Verify and integrate scattering data evaluation to the case of fully asymmetric
%   boundary conditions. 
%
%   2)  Look into implementing more efficient and stable numerical methods
%   for evaluation of scattering data. There has been intriguing research
%   into application of FFT-type methods for this purpose which,
%   apparently, maintains the symmetry between FFT and inverse FFT as between the
%   evaluation of direct scattering data to indirect scattering
%   transformation. See paper by Frumin
%   (https://www.dropbox.com/s/v419382rif7sz9t/Frumin_Efficient_Numerical_Method_for_solving_the_direct_ZS_Scattering%29_Problem_2015.pdf)
%   re. an efficient numerical algorithm for solving 
%   the inverse scattering problem and 'inverting' the method to solve the 
%   direct scattering method. Look into applying the Frumin's 'inversion' 
%   connection to apply the 'FFT'-like methods of Civelli et al
%   (https://www.dropbox.com/s/795caytlw0hkcx9/Civelli_Numerical_Methods_Inverse_Nonlinear_Fourier_Transform.pdf)
%   to the direct scattering problem for the DNLS 
%% Timing
clc;
clear all;
t = cputime;
warning('off','all'); % turns off annoying and presumably irrelevant warnings
%********************************************************************
%% DriverFile: this file ... invokes DNLS_Reader Method and passes Trial values 
    DriverFile = mfilename();
%% NumericalMethodFile: loads archived trial parameters and B-field data for analysis
    NumericalMethodFile = 'DNLS_Reader_v7'; %without '.m' extension
    DNLS_Reader_Handle = str2func(NumericalMethodFile);
%% File Locations (Automatic)
%   Call to DNLS_File_Finder() to determine
%   [ArchiveBase, ActiveArchiveBase, ParentFolder]

[ArchiveBase, ActiveArchiveBase, ParentFolder] = DNLS_File_Finder();

% ActiveArchiveBase = [ArchiveBase,'\Bradys_Trials'];
% ActiveArchiveBase = [ArchiveBase,'\Nathans_Trials'];

%NumericalMethodFile: takes trial parameters and advances solution in time
NumericalMethodFile = 'DNLS_v4'; %without '.m' extension
DNLSHandle = str2func(NumericalMethodFile);
%% System Properties: determines number of monitors at start up and uses this to open graphs on Screen #2 if available
%  If a monitor is added or removed after Matlab is started, Matlab needs to be restarted for the change to be put in effect.
            SystemGraphics = groot;
            NMonitors = length(SystemGraphics.MonitorPositions(:,1)); % returns the number of connected monitors.
            % At the moment it assumes Matlab is running on Screen #1 and,
            % if NMonitors > 1, it displays graphs on Screen #2
            
%% *************************************************************************
%% CHOOSE TRIAL # & RUN # TO ANALYZE AS WELL AS OUTPUT OPTIONS BELOW ...
%% *************************************************************************

    %% TrialNumber = designated number to which many related runs belong

    TrialNumber = '9';% The trial number stored in your Archive to analyze
    RunNumber = '9';  % ... well, the name pretty much says it all ...
    LineNumberStart = 95;   % B-field line number to study
    % LineNumberLast  = 1;    % Normally ... set to LineNumberStart
    LineNumberStepSize = 1;
    LineNumberLast = LineNumberStart; % <-- Use this as default ...
    % UNLESS evaluating real eigenvalues at a sequence of lines

    % ************************************************
    % The loop below is helpful for evaluating a sequence of real eigenvalues
    % It potentially could be helpful in automating other analysis tasks,
    % but for the moment ...TURN OFF graphing and videos and complex eigenvalue searches
    % Set EVZoom = 0 and boolREV = 1
    % ************************************************

    % When I re-ran the short Total time until coalescence, the L was set
    % to L = 15, N = 2^11 which seemed to be adequate
    ReEVidx = 0;
    XLeft = -10; % Default values
    XRight = +10; 
    for LineNumber = LineNumberStart:LineNumberStepSize:LineNumberLast 
        ReEVidx = ReEVidx + 1;

        
        %% CREATE B-FIELD GRAPH, HODOGRAM OR VIDEO?   
        boolBgraph = 1; % '1' to graph to Display B-field
        boolHodogram = 0;   % '1' to display Hodogram
        boolHodogramVideo = 0;   % '1' to display Hodogram
        boolVideo = 0; % set to 0 if boolBgraph = 1 or boolHodogram = 1
                       % set to 1 to create an avi format video
                       % set to 2 to create videos in both avi and mp4 formats
        boolIz = 0;    % used interactively to evaluate Iz across user specified section of a profile.
                       % allows user to initially zoom in and confirm range for
                       % study prior to Iz evaluation.
        boolPeakOscAnalysis = 0;
        boolSpeedCheck = 0; % designed for determining speed of a pure Fast Shock profile
                            % used in conjunction with boolVideo = 1 and will
                            % output: Vcalc, Vactual - Vcalc, Vactual
        boolDevelopment = 0;
        
        %% OPTIONS FOR EIGENVALUE SPECTRUM ANALYSIS   
            EVZoom = 1; %  = '+1' if user is to chose a smaller section of B-field for Real EV evaluation
                        %This can improve accuracy and speed up calculations
            boolEnergy = 0; % Evaluate total magnetic energy stored in profile
                            % Energy = Integral( |b|^2 - |bo|^2)dx
            boolFE = boolEnergy; % Evaluate FE = Integral{(Hilbert[|b|^2])*d/dx(|b|^2)}dx
            boolConservedQuantities = 0; % Evalautes conserved quantities for
                                         % line # = LineNumberStart to LineNumberLast
                                         % in step size = LineNumberStepSize
                            
        %% EVALUATION OF ALL REAL EIGENVALUES?        
        boolREV = 0; % set to 1 to evaluate real eigenvalues using    
            %Evaluate real eigenvalues  
                LRWin = [0.05, 0.95]; %future development: allow graphical display of region w/o evaluating entirely
                LREval = [.01, .99];  %range of real values to study
                NREval = 200;         % number of points in range to evaluate
                NBSIter = 6;          % '6' appears satisfactory and should return lr with accuracy of ...
                                      %     (lmax - lmin)/NREval/2^6 = (1 - 0)/50/64 ~ +/- 0.0003
        
        %% EVALUATION OF A SINGLE COMPLEX EIGENVALUE? 
        boolRoot = 0; %Use Newton's Method to estimate complex eigenvalue given initial guess
        zinit = 0.928+1i*0.0146; %initial guess for Newton's Method if boolRoot = 1
        
        %% EVALUATION OF ALL COMPLEX EIGENVALUES?     
        %If attempting to evaluate all Complex Eigenvalues in region LComplexEval
        boolCEV = 0; % evaluate complex eigenvalues --> Numerically intensive
        % evaluation of complex eigenvalue plane using 
        %an algorithm to identify location where 'four colors meet', these are 
        %tabulated and each is used as an initial value for Newton's method and 
        %the corresponding results are reported
                %Evaluate complex eigenvalues
                LComplexRWin = [0.01, 0.8; -.05, .99]; % future development ... as with LRWin
                LComplexEval = [0.9, 0.95; 0.0, .1]; %range of real values to study ... [Real Min, Real Max; Im Min, Im Max]
                NCEval = 20;  % Max # of intervals to evaluate along either Re(lambda) or Im(lambda) directions
        
        %% VISUAL EVALUATION OF COMPLEX EIGENVALUES?  
        
        %This option allows the user to visually verify and evaluate Complex Eigenvalues
        %and to zoom in for greater precision in evaluation
        boolCEVgraph = 0; % plot sign matrix for s11 in complex plane
     
        %% EVALUATION OF THE CONTINUOUS SPECTRUM?     
        %Necessary for representing any linear component of the wave profile
        %It has been argued that this component of the Scattering Data
        %must be responsible for, or represent, the wave profile features
        %associated with any net field rotation that is 'non-planar'
        %Additionally, if there is a change in the magnitude of the field
        %across the wave profile this must also be represented within the
        %Continuous Spectrum
        %In short: any MHD Shock must have its intrinsic features of change
        %in field magnitude or field rotation represented in the Continuous
        %Spectrum. This connection is an open question in the literature.
        
        %This is numerically intensive and needs further development
        boolRAD = 0;        % evaluate continuous scattering data
        boolRADgraph = 0;   % plot continous scattering data
        
    %% *************************************************************************
    
    %% CREATES BoolCon and EvalCon to pass User Parameters to DNLS_Reader
        boolKeys    = {'title'       , 'boolBgraph', 'boolHodogram', 'boolVideo','boolHodogramVideo', 'EVZoom', 'boolEnergy', 'boolFE', 'boolConservedQuantities', 'boolREV', 'boolRoot','zinit', 'boolCEV','boolCEVgraph', 'boolRAD', 'boolRADgraph','boolIz', 'boolSpeedCheck', 'boolDevelopment' , 'NMonitors', 'boolPeakOscAnalysis'};
        boolValues  = {'User Choices',  boolBgraph ,  boolHodogram , boolVideo  , boolHodogramVideo ,  EVZoom ,  boolEnergy ,  boolFE ,  boolConservedQuantities ,boolREV ,  boolRoot , zinit ,  boolCEV , boolCEVgraph ,  boolRAD ,  boolRADgraph , boolIz ,  boolSpeedCheck , boolDevelopment   ,  NMonitors ,  boolPeakOscAnalysis};
        BoolCon     = containers.Map(boolKeys, boolValues);
        
        EvalKeys    = {'LRWin', 'LREval', 'NREval', 'NBSIter', 'zinit', 'LComplexRWin', 'LComplexEval', 'NCEval'};
        EvalValues = { LRWin, LREval, NREval, NBSIter, zinit, LComplexRWin, LComplexEval, NCEval};
        EvalCon     = containers.Map(EvalKeys, EvalValues);
    % Verify that Re(zinit) > 0
        Proceed = 1;
        if (real(zinit)< 0) && boolRoot
            choice = questdlg('The Re(zinit) < 0 and is not a possible result.                       ', ...
            'Want to change zinit?', ...
            'Change it','Keep it','Keep it');
            % Handle response
            switch choice
    
                case 'Change it'
                    Proceed = 0;
    
                case 'Keep it'
                    Proceed = 1;
            end
        end
    %% Call to DNLS_Reader
        if Proceed == 1
            fprintf(['Begin call to: ', NumericalMethodFile,'\n']);
            [ReEVs, EVLeft_return, EVRight_return,  PStr] = DNLS_Reader_Handle(BoolCon, EvalCon, ActiveArchiveBase, TrialNumber, RunNumber, LineNumber, LineNumberStart, LineNumberLast, XLeft, XRight )    ;
            XLeft = EVLeft_return; XRight = EVRight_return;
            if (LineNumberLast > LineNumberStart) && (boolREV == 1)
                if iscellstr(ReEVs)
                    ReEVs = [0];
                end
                ReEV_list{ReEVidx} = [LineNumber,ReEVs(:)'];
            end
            
        end
        
    end
    if (LineNumberLast > LineNumberStart) && (boolREV == 1)
        Real_EV_Excel_Archive(PStr,ReEV_list,t);
    end
    
