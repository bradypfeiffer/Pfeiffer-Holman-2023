
%% Github address
% https://github.com/Dark1psoliton/DNLS_Methods
%% ToDo 6/24/2020
%   0) Create PML for Writer intended to enforce asymptotic boundary
%      https://www.dropbox.com/s/805hm9g8nbwwnfa/Paper_SpecialIssueAndre_AntoineLorinTang.pdf?dl=0
%       ... beginning on page 11 ...
%      conditions as well as eliminate cross talk others allowed by periodic
%      boundary conditions of the numerical method.
%
%   1) Create Option to Archive B-field etc. in matfile XOR txt format 
%            
%   2) Create GUI to establish user preferences etc.
%
%   3) Spin off multiple versions of functions such as IS_ShockBoundaries
%      which is used in most of the scripts as a single, stand alone script.
%
%   4) With HybridIS, improve feedback (and automate?) re. process for
%      choosing the shock position and the ramp position to optimize run time
%      for a given value of L. Communicate estimated run time before collision
%      with ramp. Or choose optimal L to provide a given required run time.
%
%   5) Address the incorrect estimate for the RHS of a generated IS that
%      arises with dispersion around 0.5 or so in the call to
%      'Generate_Intermediate_Shock'. If truncated the resulting profile is
%      either inefficient or has a discontinuity.
%% Code Updates 5/21/2023
%  Implement automatic sequence of runs where the run iteration parameter,
%  RunNumber can be used to modify one or more trial parameters.
%  The iteration loop begins on, or near, line #264 where the first
%  parameter depending on RunNumber is defined
%  The iteration loop ends on, or near 2456, after the call to
%  DNLS_Bfield_Header_v2 which is the last use of RunNumber
%  NOTE: this looping over RunNumber will not work smoothly with Profile 10
%  (a continuation of a previous run where user interaction is required) or
%  with Profile 12 (allowing the insertion of a soliton into a section of
%  an archived run ... which also requires user input).
%% Code Updates 6/29/2022
% 1) Slight update to the KBW profile. Look at comments in that section for
% guidance on selection of Iz for given compression ratio and usage in
% generating nonplanar intermediat shocks. It is hoped that this will be somewhat efficient
% for rotations near +/- 180 degrees. 
% Of course for planar intermediate shocks use the highly efficient 'Generate13IS' profile.
%
% --> changes were made to both Call_DNLS as well as DNLS_v4
% 2) The addition of the 'Continuation' profile is still in development, so don't use this yet. 
% Eventually, this will require upload of these files:
% From Writer Folder:                 Call_DNLS_v6, DNLS_v4
% From DNLS Matlab Utilities Folder: 
% Continuation_Read_Bfield_Line and Read_Bfield_Line
%
% 3) Update of 'Hybrid Profile: 1p Soliton spliced into existing Trial' ...
%   Upload to Writer Folder   : Call_DNLS_v6 & DNLS_v4
%   Upload to Utilities Folder: Continuation_Read_Bfield_Line & Read_Bfield_Line
  clear; %clear all variables from workspace
  clc;   % clear command window
%% Code Updates 5/28/2018
% 1)  added beta feature to pull Archive and Parent folder paths automatically
% from a text file created by the user and to be placed in the parent
% folder. This is intended to remove requirement on user of updating local folders
% whenever code is updated.
%
% 2) added Circular Driver where the amplitude, wavenumber, active region
% and width of transition between passive and driven regions is specified
% by user. Driver information is displayed during Writer evaluation and is
% also passed through PStr to Numerical method and, finally, these driver
% parameters are stored in the Trial Excel Archive.
%% Timing
t = cputime;
warning('off','all'); % turns off annoying and presumably irrelevant warnings
SysTime = datetime('now');
SysTime1 = rem(now,1);
% SysTime2 = rem(now,1);
% SysTimeHours = floor(SysTime2*24);
% SysTimeMinutes = floor( (SysTime2*24 - SysTimeHours)*60 );
% SysTimeSeconds = floor( ( (SysTime2*24 - SysTimeHours)*60 - SysTimeMinutes )*60 );
%********************************************************************

    
%********************************************************************  
%% File Locations (Automatic)
    Ftitle = 'File Locations';
%   Call to DNLS_File_Finder() to determine
%   [ArchiveBase, ActiveArchiveBase, ParentFolder]

    [ArchiveBase, ActiveArchiveBase, ParentFolder] = DNLS_File_Finder();

    %NumericalMethodFile: takes trial parameters and advances solution in time
    NumericalMethodFile = 'DNLS_v4'; %without '.m' extension
    DNLSHandle = str2func(NumericalMethodFile);
    
    %**********************************************************************  
%% System Properties: determines number of monitors at start up and uses this to open graphs on Screen #2 if available
            SystemGraphics = groot;
            NMonitors = length(SystemGraphics.MonitorPositions(:,1)); % returns the number of connected monitors.
            % At the moment it assumes Matlab is running on Screen #1 and,
            % if NMonitors > 1, it displays graphs on Screen #2
%% Graphics 
    Gtitle = 'Graphics Parameters';
    exact = 0;% display exact solution = 1, do not display = 0
    zmax = 1.5;
    zmin = -1.5;
    
    GraphicsKeys = {'title', 'DisplayExactSolution','zmax','zmin', 'NMonitors'};
    GraphicsValues = {Gtitle, exact, zmax,zmin, NMonitors};
    Graphics = containers.Map(GraphicsKeys, GraphicsValues);
%********************************************************************
%% !!! Numerics
    Ntitle = 'Numerics Parameters';
    N = 2^13;
%     L = 2^15*1;
    L = 10000; % 
    dx = 2*L/(N-1);
    Courant = 0.5;
    dt = dx * Courant;
    TotalTime = 200000;
    TimestepsAttempted = round(TotalTime/dt)+1;
    NumberLinesSaved = 100;
    timestepsPerFrame = round(TimestepsAttempted/(NumberLinesSaved - 1) - .5);
    
    Lva0 = 0.0;
    L2nd = 0.0;
    
    NumericsKeys   = {'title', 'N', 'L', 'dx', 'Courant', 'dt', 'timestepsPerFrame', 'TimestepsAttempted', 'TotalTime', 'NumberLinesSaved','Lva0', 'L2nd'};
    NumericsValues = {Ntitle ,  N ,  L ,  dx ,  Courant ,  dt ,  timestepsPerFrame ,  TimestepsAttempted ,  TotalTime ,  NumberLinesSaved , Lva0 ,  L2nd};
    Numerics = containers.Map(NumericsKeys, NumericsValues);

%********************************************************************

%% Physics
    Ptitle = 'Physics Parameters';
    
    %Boundary Conditions: b = by + ibz --> bo (real valued) for x --> +/- L
    bo = 0.25;
    
    %Plasma Properties
    nu = 0.05;       % dissipation
    gamma_Landau = 0.05;    % Landau damping strength
    dispersion = 0.1; % dispersion
    alpha = 0.25;      % plasma kinematics parameter
    
    %default Physics values to archive. Modified if additional physics
    %parameters required such as: 
    %   asymmetric boundary conditions, inclusion of additional physical
    %   effects (nonlinear Landau Damping, Drivers etc.)
    PhysicsKeys   = {'title','bo','nu', 'gamma_Landau', 'dispersion', 'alpha'};
    PhysicsValues = {Ptitle , bo , nu ,  gamma_Landau,   dispersion ,  alpha };
    Physics = containers.Map(PhysicsKeys, PhysicsValues);
%% !!! Choose ProfileType = 1, 2, 7, 8, ..., 12 & Trial Objective    
    ProfileType = 1; % Profiles #1 --> #12 defined below
    TrialObjective = ['Study generated 13 dispersive shock code'];
%   Stand Alone Profiles: Prfl Keys and Values Declared within profile definitions
%     ProfileType #1  = '1psoliton'; 
%     ProfileType #21 = 'double 1psoliton';
%     ProfileType #2  = '2psoliton';
%     ProfileType #3  = 'KNS';  % after 'Kaup-Newell Soliton' valid for b --> 0 at infinity
%     ProfileType #4  = 'Gaussian';
%     ProfileType #5  = 'Helical';  
%     ProfileType #6  = 'SwitchOnShock';
%     ProfileType #7  = 'FastShock';  
%     ProfileType #8  = 'KBW';  %Used by itself to form non-planar Intermediate Shocks (after Kennel, Blandford & Wu)
%     ProfileType #9  = 'GeneratedIS';  % Generates coplanar intermediate shocks

%   Used to continue an archived trial starting at a specified line number
%     ProfileType #10  = 'Continuation';

%   Hybrid Profiles: 
%   Splice a soliton into a generated shock OR into aspecified line of an archived trial
%     ProfileType #11 = Soliton spliced into a generated 13 or 23 Shock
%     ProfileType #12 = Soliton spliced into an archived profile
%% !!! Soliton Parameters: used in ProfileTypes 1, 2, 3, 11 and 12
    Soliton.Type  = 1   ;  % 1 --> 1p soliton, 21 --> double 1p soliton, 2 --> 2p soliton, 3 --> KNS soliton

    % for 1p soliton
        Soliton.lambda = 0.03;  % For a 1p soliton & actual eigenvalue is lambda_actual = lambda*bo
        Soliton.sign  = -1  ;  % +1 --> bright 1p soliton, -1 --> dark 1p soliton

    % for double 1p soliton
        Soliton.lambda1 = 0.825;  % For a 1p soliton & actual eigenvalue is lambda_actual = lambda*bo
        Soliton.sign1  = 1  ;  % +1 --> bright 1p soliton, -1 --> dark 1p soliton
        Soliton.X1  = +0.0*L  ;  % +1 --> bright 1p soliton, -1 --> dark 1p soliton
    
        Soliton.lambda2 = 0.95;  % For a 1p soliton & actual eigenvalue is lambda_actual = lambda*bo
        Soliton.sign2  = -1  ;  % +1 --> bright 1p soliton, -1 --> dark 1p soliton
        Soliton.X2  = -2.5  ;  % +1 --> bright 1p soliton, -1 --> dark 1p soliton

    % for 2p soliton
        Soliton.lr    = 0.5 ;  % Real part of eigenvalue for a 2p soliton --> actual eigenvalue is lambda_actual = (lr + I*li) *bo
        Soliton.li    = 0.5 ;  % Imaginary part of eigenvalue for a 1p soliton 

    Soliton.Locate = -1; % -1 --> place soliton on LHS of shock, +1 --> place soliton on RHS of shock
%% !!! Shock Parameters: used in ProfileTypes 7, 11 and 12
    Shock.Type  =  1;  % 1 --> 13 IS, 2 --> 23 IS on lower branch, 3 --> 23 IS on upper branch
    Shock.CR    = -0.4 ;  % Compression ratio of 13 IS with -0.5 < IS.CR13 < 0 and CR23 = -1 - IS.CR13
                           % For Fast Shock 0 < CR < 1
    Shock.Iz    = -10; % Target Iz for IS
    
    
    % theta_L and theta_kbw are only used in the KBW profile
    Shock.theta_L = 0*pi; % Normally set at 0 to match a Bfield in the vertical direction at the far LHS (by > 0 and bz = 0)
    
    Shock.theta_kbw = 180/180*pi*(- sign(Shock.Iz) ); % Direction of Bfield downstream of shock (only used in KBW)
    
    if (ProfileType == 9) || (ProfileType == 11) % Both of these profiles use GeneratedIS (which requires a planar profile)
                                                 % ... and the KBW profile to create a ramp, therefore Shock.theta_kbw must be 
                                                 % set to +/- Pi to be consistent with the planar profile to match GeneratedIS
        Shock.theta_kbw = 180/180*pi*(- sign(Shock.Iz) );
    end
    Shock.theta_R = 2.0*pi*sign(Shock.theta_kbw);
    
    % Specify location for shock and ramp
    Shock.x     = -0.35*L;  % location of midpoint of Shock in range: -L --> +L
    Shock.xRamp = +0.61*L;  % location of ramp ... in range: -L --> +L 
    
    % width of ramp on RHS of KBW profile
    Shock.rampWdth = 15;
%% !!! Trial to Create, Overwrite and Archive Trial for ProfileType 12 

    CreateMatfile = 1;  % Whether ... '1' ... or not ... '0' ... to create a matfile for archiving
                        % Note that the matfile and txt file archives are
                        % redundant and matfile usage can add significant
                        % time for trials with multiple runs. Alas
    TrialNumber = '9';  % This will create a folder, for example 'Trial 1' ... 
                        % for TrialNumber = 1 ... within folder
                        % 'Bradys_Trials' intended to store results from related trial runs
    RunNumber_start = 3;% the run number to be created with the 'Trial_N' folder
    % RunNumber_end   = 11;
    RunNumber_end   = RunNumber_start; % sets values for start/end limits on for-loop surrounding call to Numerical Method
                        % The 'RunNumber_start' and 'RunNumber_end'
                        % variables are struncture to allow multiple runs
                        % to be executed automatically where the iteration
                        % parameter called in line # 252 below can be used
                        % to modify the value of a trial parameter.
    
    % Archived B-field: Used in ProfileType 12 ... Archived Trial vs. Soliton
    TrialNumber1 = '60';     % Trial # in single quotes
    RunNumber1 = '4';        % Run # in single quotes (used in Hybrid Profile)
    LineNumber1 = 50;        % Line # as ... number ... no quotes (used in Hybrid Profile) 
%% Circular Driver Options
Drivertitle = 'Driver Parameters';

IsDriven = 0; % 1 --> turn on circular driver
DriverAmplitude = 0.002;
K_driver = 0.04;       % Driver wavenumber
X_driver = 0.25*L;   % Driver applied only in range (-X_driver, + X_driver)
W_Driver = 0.02*L;  % Width of transition between undriven & driven sections

% Default Keys & Values for Driver Parameters
DriverKeys   = {'title' ,     'IsDriven'};
DriverValues = { Drivertitle , IsDriven };
if IsDriven
    DriverKeys   = {'title' ,     'IsDriven','Amplitude',     'K',      'X',      'W'};
    DriverValues = { Drivertitle , IsDriven , DriverAmplitude, K_driver, X_driver, W_Driver };
end
Driver       = containers.Map(DriverKeys,DriverValues);
%% Nothing to see here, move along
    Ttitle = 'Trial Identifiers';
    Proceed = 0;   % Default to 0, will be set to 1 if profile construction is succesful
    BFieldarchive = 1;
    runNumericalMethod = 1;
    Excelarchive = 1;
    PlaySound = 0;  % plays sound file upon completion
    SendEmail = 0;
    CompareTrials = 0;% if == 1 --> this graphs two B-field lines allowing visual comparison
                      %         --> shifts the graph with greater N so that relevant sections coincide
                      %         --> the 'shift' feature assumes the profile
                      %         has a region defined by (EVL, EVR) ... such
                      %         as the Helical Profile
%% THE HEART OF THE CODE ... embedded with a loop over RunNumber
% with RunNumber going from RunNumber_start up to RunNumber_end
% Set RunNumber_start = RunNumber_end if you do not want to loop through
% multiple RunNumbers automatically.

for RunNumber = RunNumber_start:RunNumber_end

    %% Define dependence of a trial parameter in terms of RunNumber
        % Be sure to update any reference to run # dependent parameters
        % passed to the Numerical Method and to Archives
        if RunNumber_end > RunNumber_start
            Lva0 = 6.3;
            L2nd = -0.58;
            LMin = 5.9;
            LMax = 6.3;
            Lva0 = LMin + (RunNumber - RunNumber_start )*(LMax - LMin)/(RunNumber_end - RunNumber_start);
            NumericsKeys   = {'title', 'N', 'L', 'dx', 'Courant', 'dt', 'timestepsPerFrame', 'TimestepsAttempted', 'TotalTime', 'NumberLinesSaved','Lva0', 'L2nd'};
            NumericsValues = {Ntitle ,  N ,  L ,  dx ,  Courant ,  dt ,  timestepsPerFrame ,  TimestepsAttempted ,  TotalTime ,  NumberLinesSaved , Lva0 ,  L2nd};
            Numerics = containers.Map(NumericsKeys, NumericsValues);
        end
    %% Finalize Folder & Trial Parameters

        %Folder Paths Relative to Archive path
        TrialFolderRel = ['\Trial ',TrialNumber];%e.g. '\Trial 11'
        BfieldFolderRel = [TrialFolderRel,'\Bfields'];%e.g. '\Trial 11\Bfields'
        ArchiveBaseOld = ActiveArchiveBase;

        % Add GUI to allow for selection of reference archive, Trial # and Run #
    %     if Bradys == 1
    %         Sname = getname(Bradys);
    %         ArchiveBaseOld = ['C:\Users\rhamilto\Dropbox\Matlab Code Sandbox\Archive\',Sname,'_Trials']
    %     end
        %Absolute Folder Paths
        TrialFolderAbs = [ActiveArchiveBase,TrialFolderRel];%e.g. 'F:\ ...\Trial 11\'
        BfieldFolderAbs = [ActiveArchiveBase,BfieldFolderRel];%e.g. 'F:\ ...\Trial 11\'

        %Construct TrialID, MatID, TrialName & TimeStamp
        TrialID   = strcat(TrialNumber,'_',int2str(RunNumber)); % "#-#' format
        MatID   = strcat(TrialNumber); % "#-#' format
        TrialName = [NumericalMethodFile,'_',TrialID]; % 'DNLS2D10c_11_4'
        TimeStamp = datestr(clock,0);

        % Address for B-field *.txt file, *.mat file
        TrialK_B_field_rel = ['\Bfield_',TrialID,'.txt'];%
        B_field_Trial_Run_matfile_rel = ['\Bfield_',MatID,'.mat'];%
        B_field_Tr_Rn_txt     = [BfieldFolderAbs, TrialK_B_field_rel];
        B_field_Tr_mat     = [BfieldFolderAbs, B_field_Trial_Run_matfile_rel];
        % Open new Matfile: https://www.mathworks.com/help/matlab/ref/matfile.html
        % Matfile already exists & delete file: https://www.mathworks.com/help/matlab/ref/matfile.html
        %********************************************************************

        %Check if Trial and B-field folders exist & create if not
        utilities(TrialFolderAbs, TrialFolderRel); %Creates Trial folder if needed
        utilities(BfieldFolderAbs, BfieldFolderRel); %Creates Bfield folder if needed

        %path to excel file to record trial parameters
        xlsfilename = [NumericalMethodFile,'_',TrialNumber,'.xlsx'];% e.g. 'DNLS2D10c_11.xls'
        xlsfile = [TrialFolderAbs,'\',xlsfilename];%full path to xls file

        % Declare Trial and Folder Containers
        TrialKeys   = {'title', 'TrialNumber','RunNumber', 'TrialName','TimeStamp','TrialObjective','TrialID'};
        TrialValues = {Ttitle ,  TrialNumber , int2str(RunNumber),  TrialName,  TimeStamp , TrialObjective , TrialID};
        Trial = containers.Map(TrialKeys, TrialValues);

        FolderValues = {Ftitle,  ActiveArchiveBase ,  TrialFolderAbs,  BfieldFolderAbs,  B_field_Tr_Rn_txt,  CreateMatfile ,  B_field_Tr_mat  };
        FolderKeys = {'title' , 'ArchiveBase'      , 'TrialFolder'  , 'BfieldFolder'  , 'Bfieldfile'      , 'CreateMatfile', 'BfieldMatfile'  };
        Folder = containers.Map(FolderKeys,FolderValues);
    %********************************************************************
    %% Default Prfl Parameters  
        Prfltitle = 'Profile Parameters';
        Profile = 'Default';
        v = alpha*bo^2;   % set wave frame speed by default to linear Alfven wave speed
                    % otherwise, within a given profile, set v to soliton wave
                    % speed (or speed of slowest soliton of N-soliton
                    % profile?).

        % Default Keys & Values for Profile Parameters
        PrflKeys   = {   'title' , 'Profile', 'speed', 'bo'};
        PrflValues = { Prfltitle ,  Profile ,    v   ,  bo}; 
    %%
    %% *******************************************
    %%    Various Profile Definitions:
    %% *******************************************
    %%
    %% ProfileType #1      1p-soliton profile
    if (ProfileType == 1) || ( (ProfileType == 11 || ProfileType == 12) &&  (Soliton.Type == 1)  )

    %  if ProfileType is for a stand alone  1p-soliton ...
    %  if a 1p soliton is to be spliced into either hybrid type profile

        Profile  = '1psoliton';
        lambdaChar = char(955); % Used to print Greek letter lambda
        lambda   = Soliton.lambda;  % related to actual eigenvalue by lambda_actual = lambda*bo
        sign1 = Soliton.sign;    % sign = +1 ... 'bright soliton'; sign = -1 ... for a 'dark soliton'
        prfmsg = sprintf('Defining 1psoliton profile with %s = %4.4f and sign = %4.0f',lambdaChar,lambda, sign1);
        prfmsg

        % Prfl is a container holding the unique values required to specify properties
        % of the chosen profile. It can be passed as a single variable: Prfl
        % If these unique values are changed later, then PrflValues and Prfl
        % will need to be updated subsequent to such changes.
        % See ReadMe.docx [E:\Dropbox\Research\Kawata&Inoue\Kawata_Inoue rotationally asymmetric solitons\Dispersion is 1] for scaling transformation for case of alpha, dispersion != 1

        %*********************************************************
        v = alpha*bo^2*(1  + 2*lambda^2);  % the soliton speed ... to be updated with change 
                                % to lambda or bl & just prior to update of Prfl
        PrflKeys   = {'title', 'Profile', 'speed', 'lambda',  'sign'    , 'b1p', 'Lscale'};
        PrflValues = { Prfltitle,   Profile ,  v ,  lambda ,   sign1 ,  bo  ,  0.1*L  };
        Prfl1p       = containers.Map(PrflKeys,PrflValues);   
        %*********************************************************

        % 1p soliton definition
        c = @(Prfl, Physics) abs(Prfl('b1p'))^2 + 2*( Prfl('lambda') * Prfl('b1p') )^2;
        eta = @(Prfl, Physics) sqrt(abs(Prfl('b1p'))^2-( Prfl('lambda') * Prfl('b1p') )^2);
        OnePfactor = @(Prfl, Physics) 2*Prfl('sign')*(eta(Prfl, Physics)/(abs(Prfl('b1p'))))^2;
        Csh = @(Prfl, Physics, x,t) cosh(2*(x-c(Prfl, Physics)*t)*( Prfl('lambda') * abs( Prfl('b1p') ) )*eta(Prfl, Physics));
        Snh = @(Prfl, Physics, x,t) sinh(2*(x-c(Prfl,Physics)*t)*(( Prfl('lambda') * abs( Prfl('b1p') ) )*eta(Prfl, Physics)));
        OnePden = @(Prfl, Physics, x,t) (   Csh(Prfl, Physics, x,t)-Prfl('sign')*( Prfl('lambda') * abs( Prfl('b1p') ) )/abs(Prfl('b1p'))   ).^2;
        numy = @(Prfl, Physics, x,t) Csh(Prfl, Physics, x,t)-Prfl('sign')*abs(Prfl('b1p'))/( Prfl('lambda') * abs( Prfl('b1p') ) );
        By = @(Prfl, Physics, x, t) ( abs(Prfl('b1p')) + (( Prfl('lambda') * abs( Prfl('b1p') ) )*OnePfactor(Prfl, Physics))*(numy(Prfl, Physics, x,t))./OnePden(Prfl, Physics, x,t));
        Bz = @(Prfl, Physics, x, t) (-1)*eta(Prfl, Physics)*OnePfactor(Prfl, Physics)*(Snh(Prfl, Physics, x,t)./OnePden(Prfl, Physics, x,t));

        onePsoliton = @(Prfl, Physics, x, t) (By(Prfl, Physics, x*Physics('alpha')/Physics('dispersion'), t*Physics('alpha')^2/Physics('dispersion')) + 1i*Bz(Prfl, Physics, x*Physics('alpha')/Physics('dispersion'), t*Physics('alpha')^2/Physics('dispersion')))*(Prfl('b1p')/abs(Prfl('b1p')));

        debug1pgraph = 0;
        if debug1pgraph
            X1p = linspace(-L,L,N);
            BFieldWindowPosition = [(NMonitors - 1 +0.3) 0.65 0.4 0.3];
            BFieldFigName = sprintf('%s Graph', Prfl1p('Profile'));
            BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition,'visible','off');
            BFieldPlot(BFieldFig,onePsoliton(Prfl1p,Physics,X1p,0),X1p);
        end





        zmax = 2*1.5*max(abs(onePsoliton(Prfl1p, Physics, linspace(-20,20,500),0)));
        zmin = -zmax;
        GraphicsValues = {Gtitle, exact, zmax,zmin, NMonitors};
        Graphics = containers.Map(GraphicsKeys, GraphicsValues);

        %*********************************************************
        % Define Parameter structure & profile(Prfl, Physics, x, t) for this profile
        if (ProfileType == 1)
            Prfl = Prfl1p;
            PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
                'PrflCon', Prfl, 'GCon',Graphics,'FCon',Folder, 'DrvrCon',Driver);
            profile = @(Prfl, Physics, x, t) (onePsoliton(Prfl, Physics, x, t));
        end
        %*********************************************************

        Proceed = 1;
    end
    %******************************************************************** 
    %% ProfileType #21      double 1p-soliton profile
    if (ProfileType == 21) || ( (ProfileType == 11 || ProfileType == 12) &&  (Soliton.Type == 1)  )

    %  if ProfileType is for a double  1p-soliton ...

        Profile    = 'double_1psoliton';
        lambdaChar = char(955); % Used to print Greek letter lambda
        lambda1    = Soliton.lambda1;  % related to actual eigenvalue by lambda_actual = lambda*bo
        sign1      = Soliton.sign1;    % sign = +1 ... 'bright soliton'; sign = -1 ... for a 'dark soliton'
        X1_1p = Soliton.X1;
        beta1 = acos(lambda1/bo);

        lambda2   = Soliton.lambda2;  % related to actual eigenvalue by lambda_actual = lambda*bo
        sign2     = Soliton.sign2;    % sign = +1 ... 'bright soliton'; sign = -1 ... for a 'dark soliton'
        X2_1p     = Soliton.X2;
        prfmsg1    =   fprintf('Defining double 1psoliton profile with %s = %1.4f, sign1 = %1.0f, x1 = %4.2f\n',[lambdaChar,'1'],lambda1, sign1,X1_1p);
        prfmsg2    =   fprintf('Defining double 1psoliton profile with %s = %1.4f, sign2 = %1.0f, x2 = %4.2f\n',[lambdaChar,'2'],lambda2, sign1,X2_1p);
        beta2 = acos(lambda2/bo);

        % Prfl is a container holding the unique values required to specify properties
        % of the chosen profile. It can be passed as a single variable: Prfl
        % If these unique values are changed later, then PrflValues and Prfl
        % will need to be updated subsequent to such changes.
        % See ReadMe.docx [E:\Dropbox\Research\Kawata&Inoue\Kawata_Inoue rotationally asymmetric solitons\Dispersion is 1] for scaling transformation for case of alpha, dispersion != 1

        %*********************************************************
        v1 = alpha*bo^2*(1  + 2*lambda1^2);  % the soliton speed ... to be updated with change 
        v2 = alpha*bo^2*(1  + 2*lambda2^2);                                    % to lambda or bl & just prior to update of Prfl
        PrflKeys   = {'title'   , 'Profile', 'speed1', 'speed2','lambda1','lambda2', 'beta1', 'beta2', 'sign1', 'sign2', 'x1'  , 'x2'  , 'b1p', 'Lscale'};
        PrflValues = { Prfltitle,  Profile ,  v1     ,  v2     , lambda1 , lambda2 ,  beta1 ,  beta2 ,  sign1 ,  sign2 ,  X1_1p,  X2_1p,  bo  ,  0.1*L  };
        Prfl1p_double       = containers.Map(PrflKeys,PrflValues);   
        %*********************************************************

        % double 1p soliton definition
        
        t1 = 0; t2 = 0;
        V1   = @(Prfl, Physics) Physics('alpha')*abs(Prfl('b1p'))^2*(2 + cos(2*Prfl('beta1')));
        V2   = @(Prfl, Physics) Physics('alpha')*abs(Prfl('b1p'))^2*(2 + cos(2*Prfl('beta2')));
        k1   = @(Prfl, Physics) -Physics('alpha')*abs(Prfl('b1p'))^2*sin(2*Prfl('beta1'));
        k2   = @(Prfl, Physics) -Physics('alpha')*abs(Prfl('b1p'))^2*sin(2*Prfl('beta2'));
        Chi1 = @(Prfl, Physics, x, t) k1(Prfl, Physics).*(  x - Prfl('x1') - V1(Prfl, Physics)*(t + 0)  );
        Chi2 = @(Prfl, Physics, x, t) k2(Prfl, Physics).*(  x - Prfl('x2') - V2(Prfl, Physics)*(t + 0)  );

        EXP1  = @(Prfl, Physics, x, t) exp(-1i*Prfl('beta1')   + Chi1(Prfl, Physics,x,t));
        EXP31 = @(Prfl, Physics, x, t) exp(-3i*Prfl('beta1')   + Chi1(Prfl, Physics,x,t));
        EXP2  = @(Prfl, Physics, x, t) exp(-1i*Prfl('beta2')   + Chi2(Prfl, Physics,x,t));
        EXP32 = @(Prfl, Physics, x, t) exp(-3i*Prfl('beta2')   + Chi2(Prfl, Physics,x,t));

        Sin2b1b2 = @(Prfl, Physics) (  sin(  Prfl('beta1') - Prfl('beta2')  )/sin(  Prfl('beta1') + Prfl('beta2')  )  )^2;

        D_double = @(Prfl, Physics, x, t) 1 - Prfl('sign1').*EXP1(Prfl, Physics, x, t) - ...
            Prfl('sign2').*EXP2(Prfl, Physics, x, t) + ...
            Prfl('sign1').*Prfl('sign2').*Sin2b1b2(Prfl, Physics).*EXP1(Prfl, Physics, x, t).*EXP2(Prfl, Physics, x, t);

        B_double = @(Prfl, Physics, x, t) 1 - Prfl('sign1').*EXP31(Prfl, Physics, x, t) - ...
            Prfl('sign2').*EXP32(Prfl, Physics, x, t) + ...
            Prfl('sign1').*Prfl('sign2').*Sin2b1b2(Prfl, Physics).*EXP31(Prfl, Physics, x, t).*EXP32(Prfl, Physics, x, t);




        doubleonePsoliton = @(Prfl, Physics, x, t) abs(Prfl('b1p')).*B_double(Prfl, Physics, x, t).*conj(D_double(Prfl, Physics, x, t))./D_double(Prfl, Physics, x, t).^2;


        debug1pgraph = 0;
        if debug1pgraph
            X1p = linspace(-L,L,N);
            BFieldWindowPosition = [(NMonitors - 1 +0.3) 0.65 0.4 0.3];
            BFieldFigName = sprintf('%s Graph', Prfl1p_double('Profile'));
            BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition,'visible','off');
            BFieldPlot(BFieldFig,doubleonePsoliton(Prfl1p_double,Physics,X1p,0),X1p);
        end





        zmax = 2*1.5*max(abs(doubleonePsoliton(Prfl1p_double, Physics, linspace(-20,20,500),0)));
        zmin = -zmax;
        GraphicsValues = {Gtitle, exact, zmax,zmin, NMonitors};
        Graphics = containers.Map(GraphicsKeys, GraphicsValues);

        %*********************************************************
        % Define Parameter structure & profile(Prfl, Physics, x, t) for this profile
        if (ProfileType == 21)
            Prfl = Prfl1p_double;
            PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
                'PrflCon', Prfl, 'GCon',Graphics,'FCon',Folder, 'DrvrCon',Driver);
            profile = @(Prfl, Physics, x, t) (doubleonePsoliton(Prfl, Physics, x, t));
        end
        %*********************************************************

        Proceed = 1;
    end
    %******************************************************************** 
    %% ProfileType #2      2p-Soliton Profile
        if (ProfileType == 2) || ( (ProfileType == 11 || ProfileType == 12) &&  (Soliton.Type == 2)  )
            Profile  = '2psoliton';
            lambdaChar = char(955); % Used to print Greek letter lambda
            lr = Soliton.lr;        
            li = Soliton.li; % related to actual eigenvalue by imag(lambda) = lr*bo
            b2p = bo; % placeholder for now, it is possible that the ambient field strength for a soliton may differ from bo
            prfmsg = sprintf('Defining 2psoliton profile with %s = %4.4f + i(%4.4f)',lambdaChar,lr, li);
            prfmsg

             %*********************************************************
            v = alpha*b2p^2*( 1 + 2*(lr^2 - li^2) ); % the soliton speed ... to be updated with change 
                                        % to (lr, li) & just prior to update of Prfl

            PrflKeys   = {'title', 'Profile', 'speed', 'lr',  'li', 'b2p', 'Lscale'};
            PrflValues = { Prfltitle,   Profile ,  v ,  lr ,   li ,  b2p,   L*0.2  };
            Prfl2p       = containers.Map(PrflKeys,PrflValues); 

    %         PhysicsKeys   = {'title','bo','nu', 'gamma_Landau', 'dispersion', 'alpha'};
    %         PhysicsValues = {Ptitle , bo , nu ,  gamma_Landau,   dispersion ,  alpha };
    %         Physics = containers.Map(PhysicsKeys, PhysicsValues);          
            %*********************************************************

            rez2 = @(Prfl, Physics) ( Prfl('lr') * Prfl('b2p') )^2 - Prfl('b2p')^2 - ( Prfl('li') * Prfl('b2p') )^2;
            disc = @(Prfl, Physics) rez2(Prfl, Physics)^2 + 4*( Prfl('lr') * Prfl('b2p') )^2*( Prfl('li') * Prfl('b2p') )^2;
            sq = @(Prfl, Physics) sqrt(disc(Prfl, Physics));
            U = @(Prfl, Physics) sqrt((rez2(Prfl, Physics) + sq(Prfl, Physics))/2);
            V = @(Prfl, Physics) sqrt((- rez2(Prfl, Physics) + sq(Prfl, Physics))/2);
            Lscale = @(Prfl, Physics) 1/( ( Prfl('lr') * Prfl('b2p') )*V(Prfl, Physics) +  ( Prfl('li') * Prfl('b2p') )*V(Prfl, Physics)   );

            % LAMBDA = lambda*zeta
            CapLambda = @(Prfl, Physics)(( Prfl('lr') * Prfl('b2p') ) + 1i*( Prfl('li') * Prfl('b2p') ));  lbar = @(Prfl,Physics) (( Prfl('lr') * Prfl('b2p') ) - 1i*( Prfl('li') * Prfl('b2p') ));
            zeta = @(Prfl, Physics) (U(Prfl, Physics) + 1i*V(Prfl, Physics));      zbar = @(Prfl, Physics) (U(Prfl, Physics) - 1i*V(Prfl, Physics));
            LAMBDA = @(Prfl, Physics) CapLambda(Prfl, Physics)*zeta(Prfl, Physics);
            LBAR = @(Prfl, Physics) lbar(Prfl, Physics)*zbar(Prfl, Physics);
            k = @(Prfl, Physics) -((lbar(Prfl, Physics)*zeta(Prfl, Physics) + CapLambda(Prfl, Physics)*zbar(Prfl, Physics))/(lbar(Prfl, Physics)*zeta(Prfl, Physics) - CapLambda(Prfl, Physics)*zbar(Prfl, Physics)))^2;
            z10 = @(Prfl, Physics)  Prfl('b2p')*(CapLambda(Prfl, Physics) - zeta(Prfl, Physics))/(2*LAMBDA(Prfl, Physics)*zeta(Prfl, Physics));   z01 = @(Prfl, Physics)  Prfl('b2p')*(lbar(Prfl, Physics) + zbar(Prfl,Physics))/(2*LBAR(Prfl, Physics)*zbar(Prfl, Physics));
            z10bar = @(Prfl, Physics)  Prfl('b2p')*(lbar(Prfl, Physics) - zbar(Prfl, Physics))/(2*LBAR(Prfl, Physics)*zbar(Prfl, Physics));   z01bar = @(Prfl, Physics)  Prfl('b2p')*(CapLambda(Prfl, Physics) + zeta(Prfl, Physics))/(2*LAMBDA(Prfl, Physics)*zeta(Prfl, Physics));
            d10 = @(Prfl, Physics)  ((CapLambda(Prfl, Physics) - zeta(Prfl, Physics))/Prfl('b2p'))^2 * z10(Prfl, Physics);   d01 = @(Prfl, Physics)  ((lbar(Prfl, Physics) + zbar(Prfl, Physics))/Prfl('b2p'))^2 * z01(Prfl, Physics);

            phi = @(Prfl, Physics, x,t) 2*LAMBDA(Prfl, Physics)*(x - (Prfl('b2p')^2 + 2*(( Prfl('lr') * Prfl('b2p') )^2 - ( Prfl('li') * Prfl('b2p') )^2 + 2i*( Prfl('lr') * Prfl('b2p') )*( Prfl('li') * Prfl('b2p') )))*t);
            phibar = @(Prfl, Physics, x,t) 2*LBAR(Prfl, Physics)*(x - (Prfl('b2p')^2 + 2*(( Prfl('lr') * Prfl('b2p') )^2 - ( Prfl('li') * Prfl('b2p') )^2 - 2i*( Prfl('lr') * Prfl('b2p') )*( Prfl('li') * Prfl('b2p') )))*t);

            f = @(Prfl, Physics, x, t) 1 + exp(-1i*phi(Prfl, Physics,x,t))/(Prfl('b2p')*z10(Prfl, Physics));
            g = @(Prfl, Physics, x, t) 1 + exp(1i*phibar(Prfl, Physics,x,t))/(Prfl('b2p')*z01(Prfl, Physics));

            fbar = @(Prfl, Physics,x,t) 1 + exp(1i*phibar(Prfl, Physics,x,t))/(Prfl('b2p')*z10bar(Prfl, Physics));
            gbar = @(Prfl, Physics,x,t) 1 + exp(-1i*phi(Prfl, Physics,x,t))/(Prfl('b2p')*z01bar(Prfl, Physics));

            m = @(Prfl, Physics,x,t) 1 + exp(-1i*phi(Prfl, Physics,x,t))/(Prfl('b2p')*d10(Prfl, Physics));
            n = @(Prfl, Physics,x,t) 1 + exp(1i*phibar(Prfl, Physics,x,t))/(Prfl('b2p')*d01(Prfl, Physics));

            %valid if alpha = 1 and dispersion = 1
            twopSolitonUnscaled = @(Prfl, Physics,x,t) Prfl('b2p')*(m(Prfl, Physics,x,t).*n(Prfl, Physics,x,t) - (1 + k(Prfl, Physics))).*(fbar(Prfl, Physics,x,t).*gbar(Prfl, Physics,x,t) - (1 + k(Prfl, Physics)))./(f(Prfl, Physics,x,t).*g(Prfl, Physics,x,t) - (1 + k(Prfl, Physics))).^2;

            %valid for any real, positive alpha & dispersion values
            twopSoliton = @(Prfl, Physics,x,t)  twopSolitonUnscaled(Prfl, Physics,x*Physics('alpha')/Physics('dispersion'),t*Physics('alpha')^2/Physics('dispersion'));

            zmax = 2*1.5*max(abs(twopSoliton(Prfl2p, Physics, linspace(-20,20,500),0)));
            zmin = -zmax;
            GraphicsValues = {Gtitle, exact, zmax,zmin, NMonitors};
            Graphics = containers.Map(GraphicsKeys, GraphicsValues);
            %*********************************************************
            % Define Parameter structure & profile(Prfl, Physics, x, t) for this profile
            if (ProfileType == 2)
                Prfl = Prfl2p;
                PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
                    'PrflCon', Prfl, 'GCon',Graphics,'FCon',Folder, 'DrvrCon',Driver);
                profile = @(Prfl, Physics,x,t) twopSoliton(Prfl, Physics,x,t);;
            end
            %*********************************************************
            Proceed = 1;
        end
    %% ProfileType #3      KNS
    % The 2p-soliton discovered by Kaup & Newell [J. Math. Phys. 19(4), 1978]
    % valid for boundary conditions with b --> 0 to the left and right
    if (ProfileType == 3) && 0
    %     lr  = 1*cos(5*pi/180);      % Default value for Re(lambda)
    %     li    = 1*sin(5*pi/180);    % Default value for Im(lambda)
        MagLambda = 1;  % magnitude of the complex eigenvalue
        KN_angle = 0.05; % pick 0 < KN_angle < 1 ... angle measured clockwise from the '-' axis
                        % 0.05 --> B peak of ~ 0.1
                        % 0.5  --> B peak of ~ 1


        % Prfl is a container holding the unique values required to specify properties
        % of the chosen profile. It can be passed as a single variable: Prfl
        % If these unique values are changed later, then PrflValues and Prfl
        % will need to be updated subsequent to such changes.

        %*********************************************************
        v = 4*MagLambda*cos( pi*KN_angle ); % the KN soliton speed so that xo = v*t in KN notation
                                % will need to be updated with changes to lr & li
        PrflKeys   = {'title', 'Profile', 'speed', 'MagLambda',  'KN_angle'};
        PrflValues = { Prfltitle,   Profile ,  v ,  MagLambda ,   KN_angle };
        Prfl       = containers.Map(PrflKeys,PrflValues);

    %     PhysicsKeys   = {'title','bo','nu', 'dispersion', 'alpha'};
    %     PhysicsValues = {Ptitle , bo , nu ,  dispersion ,  alpha };
    %     Physics = containers.Map(PhysicsKeys, PhysicsValues);   
        %*********************************************************

        % KN soliton definition
    %     KN_del = @(Prfl) (Prfl('lr')^2 + Prfl('li')^2)^(1/4);
        KN_del = @(Prfl) sqrt( Prfl('MagLambda') );
        gamma = @(Prfl) pi*Prfl('KN_angle');
        eta = @(Prfl) KN_del( Prfl )^2 * sin( gamma(Prfl) );
        xi = @(Prfl) KN_del( Prfl )^2 * cos( gamma(Prfl) );
    %     gamma = @(Prfl) pi - atan(Prfl('li')/Prfl('lr')); %KN define gamma from -x axis. Ya, go figure.
        xo_dot =  @(Prfl) +4*KN_del(Prfl)^2*cos(gamma(Prfl));
        so_dot = @(Prfl) 2*KN_del(Prfl)^4*cos(2*gamma(Prfl));

        theta = @(Prfl, x, t) eta(Prfl)*( x + xo_dot(Prfl)*t);
        sigma = @(Prfl, x, t) xi(Prfl)*x - so_dot(Prfl)*t;
        exppg = @(Prfl, x, t) exp( 4*theta(Prfl, x, t) ) + ( cos( gamma(Prfl) ) + 1i*sin( gamma(Prfl) ) );
        expmg = @(Prfl, x, t) exp( 4*theta(Prfl, x, t) ) + exp(-1i* gamma(Prfl) );
        num = @(Prfl, x, t) exp( 2*theta(Prfl, x, t) ).* exp(- 2i*sigma(Prfl,x,t) ) .* expmg( Prfl, x, t);
        OnePden = @(Prfl, x, t) exppg(Prfl,x,t).^2;

        KNS = @(Prfl, Physics, x, t) 4 * KN_del(Prfl) * sin( gamma(Prfl) ) .* num( Prfl, x, t ) ./ OnePden( Prfl, x, t);
        profile = @(Prfl, Physics, x, t) (KNS(Prfl, x, - t));
    end
    %********************************************************************
    %% ProfileType #4      Gaussian Profile            
    if (ProfileType == 5) && 0

        amplitude = 1;
        W = 5;
        K = 0;

        PhysicsKeys   = {'title','bo','nu', 'dispersion', 'alpha'};
        PhysicsValues = {Ptitle , bo , nu ,  dispersion ,  alpha };
        Physics = containers.Map(PhysicsKeys, PhysicsValues);   

    %     v = bo^2;   % set by default to Alfven wave frame above
    %    % W_AS = 0.6*L; %width of helix: proportional to # of solitons in profile
    %     W_AS = 50*pi; %width of helix: proportional to # of solitons in profile
    %     W_0 = W_AS; % transition to periodic BCs
    %     L_AS = -0.3*L; %center of helix
    % %     K_Helix = - 2*pi/W_AS; %helix wave number: key to identify stability regime of 4 possible cases
    %     K_Helix =  -1*pi/W_AS; %helix wave number: key to identify stability regime of 4 possible cases
    %     helix_angle = K_Helix*W_AS;   
    %     theta_p = mod(helix_angle,2*pi);
    %     xL = L_AS - (W_AS/2); EVL = (xL + (-L))*0.5; %1/2 way between xL and (-L)
    %     xR = L_AS + (W_AS/2); EVR = xR + 0.25*(L - xR);%1/4 of the way from xR to L

    %     modangle = mod(helix_angle,2*pi); %returns diff. between helix-angle and 2npi as 0 < modangle < 2*pi
    %     if modangle > pi
    %         modangle = 2*pi - modangle; % returns diff. as -pi < modangle < +pi ... so rotation is smallest possible
    %     end
    %     
    %     theta_0 = @(x,t) (1+tanh((x-(xR+0.75*(L - xR)))/W_0))*(-.5)*mod(helix_angle,2*pi);
    %     theta_0 = @(x,t) 0;
    %     theta_h = @(x,t)  (x>=xL & x<=xR).*((x-xL)/W_AS)*helix_angle + ...
    %         (x>xR).*helix_angle;   
    %     theta_h = @(x,t) (1+tanh((x-xL)/W_AS))*(.5)*helix_angle;
        By = @(x,t,amplitude, W, K) amplitude*exp(-x.*x/W).*cos(K*x);
        Bz =  @(x,t,amplitude, W, K) amplitude*exp(-x.*x/W).*sin(K*x);

        PrflKeys = {'title', 'Profile','speed', 'amplitude', 'width', 'K'};
        PrflValues = {Prfltitle, Profile , v     , amplitude ,  W   , K};
        profile = @(Prfl, Physics, x,t) (By(x,t,amplitude, W, K) + 1i*Bz(x,t,amplitude, W, K));
    end
    %% ProfileType #5      Helical Profile            
    if (ProfileType == 5) && 0

        PhysicsKeys   = {'title','bo','nu', 'dispersion', 'alpha'};
        PhysicsValues = {Ptitle , bo , nu ,  dispersion ,  alpha };
        Physics = containers.Map(PhysicsKeys, PhysicsValues);   

        v = bo^2;   % set by default to Alfven wave frame above
       % W_AS = 0.6*L; %width of helix: proportional to # of solitons in profile
        W_AS = 50*pi; %width of helix: proportional to # of solitons in profile
        W_0 = W_AS; % transition to periodic BCs
        L_AS = -0.3*L; %center of helix
    %     K_Helix = - 2*pi/W_AS; %helix wave number: key to identify stability regime of 4 possible cases
        K_Helix =  -1*pi/W_AS; %helix wave number: key to identify stability regime of 4 possible cases
        helix_angle = K_Helix*W_AS;   
        theta_p = mod(helix_angle,2*pi);
        xL = L_AS - (W_AS/2); EVL = (xL + (-L))*0.5; %1/2 way between xL and (-L)
        xR = L_AS + (W_AS/2); EVR = xR + 0.25*(L - xR);%1/4 of the way from xR to L

        modangle = mod(helix_angle,2*pi); %returns diff. between helix-angle and 2npi as 0 < modangle < 2*pi
        if modangle > pi
            modangle = 2*pi - modangle; % returns diff. as -pi < modangle < +pi ... so rotation is smallest possible
        end

        theta_0 = @(x,t) (1+tanh((x-(xR+0.75*(L - xR)))/W_0))*(-.5)*mod(helix_angle,2*pi);
    %     theta_0 = @(x,t) 0;
    %     theta_h = @(x,t)  (x>=xL & x<=xR).*((x-xL)/W_AS)*helix_angle + ...
    %         (x>xR).*helix_angle;   
        theta_h = @(x,t) (1+tanh((x-xL)/W_AS))*(.5)*helix_angle;
        By = @(x,t) bo*cos(theta_h(x,t)+theta_0(x,t));
        Bz = @(x,t) bo*sin(theta_h(x,t)+theta_0(x,t));

        PrflKeys = {'title', 'Profile','speed', 'W_AS', 'W_0', 'L_AS', 'K_Helix', 'helix_angle', 'theta_p', 'xL', 'xR','EVL', 'EVR'};
        PrflValues = {Prfltitle, Profile , v     , W_AS ,  W_0 ,  L_AS ,  K_Helix ,  helix_angle ,   theta_p ,  xL ,  xR , EVL ,  EVR};
        profile = @(Prfl, Physics, x,t) (By(x,t) + 1i*Bz(x,t));
    end
    %% ProfileType #6      Switch-On Shock
    if (ProfileType == 6) && 0

        % Prfl is a container holding the unique values required to specify properties
        % of the chosen profile. It can be passed as a single variable: Prfl
        % If these unique values are changed later, then PrflValues and Prfl
        % will need to be updated subsequent to such changes.

        %*********************************************************
        b1 = bo*.9;
        v = alpha*b1^2;  % SO-shock speed, update Prfl upon change
        SO_sign = 1;
        SO_ramp = -SO_sign;
        XSO = -L/2;
        XR = +L/2;
        WR = L/10;
        PrflKeys   = {'title'   , 'Profile', 'speed', 'XSO', 'XR', 'BSO', 'WR'};
        PrflValues = { Prfltitle,  Profile ,  v     ,  XSO ,  XR ,  b1  ,  WR};
        Prfl       = containers.Map(PrflKeys,PrflValues);

        PhysicsKeys   = {'title','bo','nu', 'dispersion', 'alpha'};
        PhysicsValues = {Ptitle , bo , nu ,  dispersion ,  alpha };
        Physics = containers.Map(PhysicsKeys, PhysicsValues);   
        %*********************************************************

        Lbar =@(Prfl) ( Prfl('speed') )./( abs(Prfl('speed') ) )*( dispersion^2 + nu^2 )/( alpha * Prfl('BSO')^2 * nu );

        transition = @(Prfl,x) 1/2 * ( 1 + tanh( (x - Prfl('XR'))/Prfl('WR') ) );

        f = @(Prfl, x, t) (1 + exp( 2 *( x - Prfl('XSO') - Prfl('speed') * t )/(Lbar(Prfl))) );
    %     fr = @(Prfl, x, t) (1 + exp( 2 *( x - Prfl('XR') - Prfl('speed') * t )/(*Lbar(Prfl))) );

        SwitchOnShock = @(Prfl, x, t) Prfl('BSO') * ( 1 ./ sqrt( f( Prfl, x, t ) ) ) .* ( cos( dispersion ./ (2* nu) .* log( f( Prfl, x, t ) ) ) - 1i* sin( dispersion ./ (2* nu) .* log( f( Prfl, x, t ) ) ) );  
    %     SwitchOnShock_ramp = @(Prfl, x, t) Prfl('BSO') * ( 1 ./ sqrt( fr( Prfl, x, t ) ) ) .* ( cos( dispersion ./ (2* nu) .* log( fr( Prfl, x, t ) ) ) - 1i* sin( dispersion ./ (2* nu) .* log( fr( Prfl, x, t ) ) ) );    

        profile = @(Prfl, Physics, x, t) SwitchOnShock(Prfl, x, t) + SwitchOnShock(Prfl, -x + Prfl('XSO') + Prfl('XR'), -t);
    end
    %********************************************************************
    %% ProfileType #7      Fast Shock Profile
    if (ProfileType == 7)

        Profile  = 'FastShock';
         %*********************************************************
            bl = bo;
            br = abs( Shock.CR )*bl;
            XFS = Shock.x;
            XRHS = Shock.xRamp;
            WFS = L/50;
    %         WFS = (nu^2 + dispersion^2)/(alpha*nu*bl^2)/(0.33152);
            WRHS = L/50;
            v = alpha*(bl^2 + br^2 + bl*br);
    %         v = alpha*(1.1563*bl^2 +1.55784*br^2+0.3834*br*b1); % the shock speed ... to be updated with change 
                                        % to (lr, li) & just prior to update of Prfl
            PrflKeys   = {'title'   , 'Profile', 'speed', 'bl',  'br', 'xfs', 'xrhs', 'wfs', 'wrhs'};
            PrflValues = { Prfltitle,  Profile ,  v     ,  bl ,   br,   XFS ,  XRHS ,  WFS ,  WRHS };
            PrflFS       = containers.Map(PrflKeys,PrflValues);
         %*********************************************************
         % Tanh Profile ... to evolve into steady Fast Shock ... eventually
         FSprofile = @(Prfl, Physics,x,t) (Prfl('bl') - Prfl('br'))*(1/2)*(1 - tanh( (x - Prfl('xfs'))/Prfl('wfs') )) + Prfl('br') -...
                              (Prfl('br') - Prfl('bl'))*(1/2)*(1 + tanh( (x - Prfl('xrhs'))/Prfl('wrhs') ));

         %*********************************************************
         % Update zmax and Graphics container
         zmax = 1.5*max(abs(FSprofile(PrflFS, Physics, linspace(-L,L,N),0)));
         zmin = -zmax;
         GraphicsValues = {Gtitle, exact, zmax,zmin, NMonitors};
         Graphics = containers.Map(GraphicsKeys, GraphicsValues);
         % Define Parameter structure & profile(Prfl, Physics, x, t) for this profile
         if (ProfileType == 7)
             Prfl = PrflFS;
             PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
                 'PrflCon', Prfl, 'GCon',Graphics,'FCon',Folder, 'DrvrCon',Driver);
             profile = @(Prfl, Physics,x,t) ( FSprofile(Prfl, Physics, x, t) );
         end
         %*********************************************************
        Proceed = 1;



    end
    %% ProfileType #8      KBW Profile
    if (ProfileType == 8) || (ProfileType == 9) || (ProfileType == 11)
        %% The formula for Wkbw produces an Intermediate shock with an Iz that
        %matches value specified on a preceeding line. If Iz is the integral of
        %Im(b) from 'minus infinity to plus infinity' ... that is across the
        %range used to model the shock development ... theta_kbw = +/- Pi, that
        %is for a co-planar shock, the integral is well defined and unique. For
        %a non-coplanar shock Im(b) does not go to zero to the right. Because
        %of this, any shock movement would change the distance integrated on
        %the downstream side of the shock which would adjust the numerical
        %result of integration of Im(b) making it's meaning problematic. Presumably
        %this might be accounted for if it is possible to uniquely identify,
        %algorithmically, the beginning and ending of the shock.
        %
        % An additional observation is that Iz, assuming there exists a way to make a well defined
        % integral of Im(b) across the shock, will not be constant in time for
        % a non-coplanar shock but must increase monotonically in time. This is
        % a critical feature of Intermediate Shock dynamics described by KBW.
        %
        % It appears possible that there <may> be another conserved quantity
        % being the integral of Re(b)*Cos(phi) + Im(b)*Sin(phi) for a specific
        % value of phi determined by the compression ratio and the downstream
        % shock angle, theta_kbw. Such an angle has been tentatively shown to
        % exist, at least in some cases. If so, it is unclear what relation it
        % has to shock dynamics.

        %% KBW Parameters
        Profile = 'KBW';
        boolShockFit = 1;

        b_L = bo;        %Left hand boundary conditions
        theta_L = Shock.theta_L;    


        b_kbw =abs ( Shock.CR )*b_L;           % Boundary conditions to right of KBW profile

        Iz = Shock.Iz;  % Look at critical curve for given dispersion, dissipation and compression ratio. 
                        % While the critical curve graphs are only valid for planar
                        % profiles, for angles close to +/- 180 they should be close
    %     theta_kbw = 180/180*pi*(-sign(Iz));

        theta_kbw = Shock.theta_kbw;% Place '-1' in parentheses for to make Iz > 0  and '+1' to make Iz < 0.

    %     Wkbw = 1.56;
        % There is a well defined connection between Wkbw and Iz for the KBW
        % profile: Given the value of Iz above, the choice of Wkbw evaluated
        % below will produce a profile which will have the specified Iz. This
        % is valid for planar profiles and should be a good approximation for
        % rotations near +/- 180 degrees.
        Wkbw = 2*abs(Iz)/((b_kbw + b_L)*abs(Si(theta_kbw))); % Evaluated using a Pade approximation as defined above. Of dubious meaning for non-coplanar profiles.
        xkbw = Shock.x;

        b_R = b_L;      %Right hand boundary
        theta_R = Shock.theta_R;
        xRkbw = Shock.xRamp;  % location of ramp
        WR = Shock.rampWdth;      % width of ramp

        %*********************************************************
        v = (1 + (Shock.CR)^2 + (Shock.CR) )*alpha*abs(b_L)^2;  % the IS shock speed                          
    %     Iz = 0; %placeholder until evaluation at final definition of PStr
        CR = -1/2*(1 + sqrt(1 + 4*( v - b_L^2))); 
        lambda3 = v/(alpha*b_L^2) - 1;
        r3 = sqrt(lambda3 + 0.25);
        q13 = -1/2 + r3;
        q23 = -1/2 - r3;
        PrflKeys   = {'title',    'Profile', 'speed', 'CR13', 'CR23', 'Wk'  , 'Wr' , 'theta_L', 'theta_kbw', 'Iz', 'theta_R' , 'b_L', 'b_kbw', 'b_R','Xkbw','Xr'   , 'boolShockFit'};
        PrflValues = { Prfltitle,  Profile ,  v     ,  q13     q23  ,  Wkbw ,  WR  ,  theta_L ,  theta_kbw ,  Iz ,  theta_R  ,  b_L ,  b_kbw ,  b_R , xkbw , xRkbw ,  boolShockFit };
        PrflKBW       = containers.Map(PrflKeys,PrflValues);
        %*********************************************************

        transition = @(x,y,w) 1/2 * ( 1 +  tanh( (x - y)/w ) );

        bT = @(Prfl,x) Prfl('b_L') + ( Prfl('b_kbw') - Prfl('b_L') ) * transition(x, Prfl('Xkbw'), Prfl('Wk')) + ...
            + (Prfl('b_R') - Prfl('b_kbw') )*transition(x, Prfl('Xr'), Prfl('Wr'));
        theta = @(Prfl,x) Prfl('theta_L') + (Prfl('theta_kbw') - Prfl('theta_L'))*transition(x, Prfl('Xkbw'), Prfl('Wk')) + ...
            + (Prfl('theta_R') - Prfl('theta_kbw'))*transition(x, Prfl('Xr'), Prfl('Wr'));

        bkbw = @(Prfl,x) bT(Prfl,x) .* ( cos( theta(Prfl,x) ) + 1i*sin( theta(Prfl,x) ) ) ;

        KBWprofile = @(Prfl, Physics, x,t)  bkbw(Prfl,x);

        zmax = 1.5*max(abs(KBWprofile(PrflKBW, Physics, linspace(-L,L,N),0)));
        zmin = -zmax;
        GraphicsValues = {Gtitle, exact, zmax,zmin, NMonitors};
        Graphics = containers.Map(GraphicsKeys, GraphicsValues);
    %%  Evaluate Speeds, Distances, Time to BC Failure
    %   The ramp on the RHS will form into an intermediate shock travelling
    %   slower than the intended shock on the LHS. 
    %   The Iz of the ramp correlates to the compression ratio of its resulting
    %   IS and, thus, the decrease in the B-field magnitude on the RHS of the
    %   ramp IS. The wider the ramp is, the bigger its Iz will be making its
    %   compression ratio closer to -1, which cause it to travel at low speed
    %   in the lab frame. The smaller the ramp's Iz, the smaller the magnitude
    %   of its compression ratio, the faster it will travel and the longer
    %   before the intended IS on the LHS will collide with it.
    %   The faster the ramp IS travels, the smaller the B-field on its RHS ...
    %   and the faster the rarefactive wave on the far RHS will travel, thus,
    %   causing it to collide with the intended IS on the LHS faster. So, there
    %   must be some optimal width, or Iz, for the ramp that will allow the
    %   trial to run for a longer time for a given value of L.
        Lramp = log(1000)*WR;
        Lshock = log(1000)*Wkbw; % should be roughly valid ... without dispersion!

        XShockR = Shock.x + Lshock/2;
        XShockL = Shock.x - Lshock/2;

        LPrfl = Lshock + Lramp; % width of shock and ramp ... 

        %Left and Right edges of Ramp
        XrampL = Shock.xRamp - Lramp/2;
        XrampR = Shock.xRamp + Lramp/2;

        % distance from left edge of ramp to right edge of 'Shock-Soliton' profile
        DXRight = XrampL - XShockR; % from RHS of shock to LHS of ramp

        % distance from right edge of ramp to RHS = L - XrampL
        % distance from left edge of 'shock-soliton' profile = min(XLsoliton, XShockL) - ( -L)
        DXLeft = ( L - XrampR ) + ( XShockL - (-L) ); % distance from RHS of ramp across boundary to LHS of shock

        DXfree = 2*L - LPrfl;

        VRWFoot = alpha*(b_kbw^2); % roughly valid if Iz of ramp is large enough
        VRW = alpha*(bo^2 + b_kbw^2 + bo*b_kbw); % Seems to travel at speed of FS

        PrfWidth = LPrfl;
        CollisionDistance = (2*L - PrfWidth  )/2;

        Timpact_Left = DXLeft/(VRW - v);
        Timpact_Right = DXRight/(v - VRWFoot);
        estImpactTime = min(Timpact_Left,Timpact_Right);

        Lnew = (LPrfl + (VRW - VRWFoot)*TotalTime)/2;
        DL_new = (VRW - v)*TotalTime;
        Xshock_new =( -Lnew + (Lshock + DL_new)/2 )/Lnew;
        Xramp_new = (Lnew - (DL_new + Lramp)/2 )/Lnew;
        NnewTmp = 2*Lnew/0.015;
        powbase2 = floor(log(NnewTmp)/log(2)) + 1;
        Nnew = 2^powbase2;
        Proceed = 0;
        ProceedKBW = 0;
        msg1 = ['Suggest the following changes ...',newline,newline];
        msg1 = [msg1,'Lnew = ',num2str(Lnew),newline];
        msg1 = [msg1,'Nnew = ',num2str(Nnew),' ... or 2^p with p = ',num2str(powbase2),newline];
        msg1 = [msg1,'Shock.x = (', num2str(Xshock_new),')*Lnew',newline];
        msg1 = [msg1,'Shock.xramp = (', num2str(Xramp_new),')*Lnew']
        if estImpactTime < TotalTime
            if estImpactTime < 0
                msg1=['BC''s fail before trial even starts :-(',newline];
            else
                msg1 = ['BC''s will fail at about t = ',num2str(estImpactTime),newline];
            end
            msg1 = [msg1,'Suggest the following changes ...',newline,newline];
            msg1 = [msg1,'Lnew = ',num2str(Lnew),newline];
            msg1 = [msg1,'Nnew = ',num2str(Nnew),' ... or 2^p with p = ',num2str(powbase2),newline];
            msg1 = [msg1,'Shock.x = (', num2str(Xshock_new),')*Lnew',newline];
            msg1 = [msg1,'Shock.xramp = (', num2str(Xramp_new),')*Lnew']

            choice = questdlg('L is too small and BC''s may fail too soon:', ...
                'Keep doomed L as is or restart and Adjust values?',...
                'Keep','Adjust', 'No');
            % Handle response
            switch choice
                case 'Keep'
                    disp([choice ' L as is. What could go wrong?'])
                    Proceed = 1; ProceedKBW = 1;
                case 'Adjust'
                    disp([choice ' with suggestions above ...'])
                    Proceed = 0; ProceedKBW = 0;
                case 'No'
                    disp([choice ' I''m out of here!'])
                    Proceed = 0; ProceedKBW = 0;
            end
        else
            Proceed = 1;
            ProceedKBW = 1;

        end


    %     [XrampNew, XISNew, TBCFail, Tcollision] = HybridOptimize(Lnew, LPrfl, Lramp, Lsoliton, LShock, Soliton, vsol, vIS, VRW, VRWFoot);

    %%     

        if ProfileType == 8 && Proceed == 1
            Prfl = PrflKBW;
            PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
                    'PrflCon', Prfl, 'GCon',Graphics,'FCon',Folder, 'DrvrCon',Driver);
            profile = @(Prfl, Physics, x,t) KBWprofile(Prfl, Physics,x,t);    
        elseif Proceed == 1
            remove(PrflKBW,'Profile');
            remove(PrflKBW,'speed');
            remove(PrflKBW,'CR13');
            remove(PrflKBW,'CR23');
            remove(PrflKBW,'Iz');
            PrflKBW('Profile3') = 'KBW';
        end

    end
    %% ProfileType #9      GeneratedIS
     if ((ProfileType == 9) || (ProfileType == 11)) &&  ProceedKBW == 1
        Profile = 'GeneratedIS';
        %% Planar Intermediate Shock Parameters and Options Specified above through 'Shock' object

        CR = Shock.CR; % Compression ratio for Intermediate Shock
                       % -1 < CR < 0   
        Iztarget = Shock.Iz;  % Iz is the area under Im(b) ... integrated from right to left
                         % So if Im(b) is positive throughout the shock, then
                         % Iz will be negative
        xIS = Shock.x;  % location of midpoint of IS Shock in range: -L --> +L
        boolIS = Shock.Type; % set to '1' to search for 1->3 Shock with Iz close to I13target
                    % set to '2' to search for 2->3 Shock with smallest Iz
                    % set to '3' to search for 2->3 Shock with largest Iz

        boolShockFit = 0;

        boolCallShockGraphing = 1; % set to 1 for debugging 'Call_DNLSB_Given_Compression_Ratio'

        boolDebugProfileGraphs = 1; % set to 1 to generate graphs to help debug splicing 1->3 Shock onto background

        b_L = bo;        %Left hand boundary conditions
        theta_L = 0*pi;

        xRkbw = Shock.xRamp;  % location of ramp
        WR = L/40;      % width of ramp

        theta_IS = 180/180*pi*(-sign(Iztarget)); % Generated13IS is currently only defined for theta_IS = +(-)Pi for Iztarget negative (positive)
        bLasymptote = 0; % This is the value of the z-component of the magnetic field far to the LEFT of the shock.
                         % Throughout (and consistent with KBW) this is taken
                         % to be zero
        bRasymptote = 0; % This is the value of the z-component of the magnetic field far to the RIGHT of the shock.
                         % For a planar profile this will be 0.
                         % As the profile generation method looks for
                         % stationary  field profiles connecting the stationary
                         % points and as non-coplanar profiles are not
                         % stationary, it is not obvious that an ODE based
                         % method can quickly generate a shock solution of
                         % non-coplanar profiles. It is possible to reduce the
                         % DNLSB to an ODE through a similarity transformation.

        b_R = b_L;      %Right hand boundary conditions
        theta_R = 2.0*pi*sign(theta_IS);

        b_kbw = abs(CR)*bo;   

        NTEST = 10; % (suggested to use 10) The number of points to average over in the search for shock boundaries

        TOL = 0.00002; % (suggested to use 0.001) The allowed discrepancy between Shock Bz and the asymptotic values defined above   
        %% Call to Generate 1->3 Shock

        if boolIS == 1
            prfmsg = sprintf('Looking for 1->3 shock profile with CR_{13} = %4.4f, Iz =  %4.4f',CR,Iztarget);
        elseif boolIS ==2
            prfmsg = sprintf('Looking for 2->3 shock profile with CR_{23} = %4.4f with the lowest Iz',CR);
        elseif boolIS ==3
            prfmsg = sprintf('Looking for 2->3 shock profile with CR_{23} = %4.4f with the highest Iz',CR);
        end 
        prfmsg

        [I13best, B13best, X13best, betaBest, ShockFound, msgBG] = Generate_Intermediate_Shocks(CR, Iztarget, bLasymptote, bRasymptote, NTEST, TOL, dispersion, nu, alpha, bo, boolCallShockGraphing, boolIS);

        if boolIS == 1 && ShockFound
            prfmsg1 = sprintf('Found 1->3 shock profile for CR_{13} = %4.4f with Iz =  %4.4f',CR,I13best);
            pctError = abs(I13best-Iztarget)/(abs(Iztarget) + 0.001)*100;
            prfmsg2 = sprintf('This Iz of %4.4f is %4.0f percent off from the target Iz of %4.4f',I13best,pctError,Iztarget);
            Iztarget = I13best;
            prfmsg3 = sprintf('Setting Iztarget = IzFound = %4.4f',I13best);
            prfmsg3
            prfmsg1
            prfmsg2

            if pctError > 10
                 % Construct a questdlg with three options
                    choice = questdlg('Proceed even though % error is large?', ...
                    'Proceed?:', ...
                    'Yes','Do not proceed','Do not proceed');
                    % Handle response
                    switch choice
                        case 'Yes'
                            disp([choice ' ... Proceeding.'])
                            ShockFound = 1;
                        case 'Do not proceed'
                            disp([choice ' Stopping now.'])
                            ShockFound = 0;
                        case 'Stopping now.'
                            disp([choice ' Stopping now.'])
                            ShockFound = 0;
                    end   
            else 
                Iztarget = I13best;
            end
        end 
        %BFieldWindowPosition = [(NMonitors - 1 +0.3) 0.65 0.4 0.3];BFieldFigName = 'Magnetic Field Graph';BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition);plot(X13best,imag(B13best))
        %% Fit 1->3 Shock Uniform Grid & update N if necessary
        if ShockFound
            msgBG
            if X13best(end) < X13best(1)
                X13best = flip(X13best);
                B13best = flip(B13best);
            end
            % Run through Shock Boundary Location function to reduce to actual
            % length of shock.
             tol = 0.001;
             ntest = 10;%(bb,x,bL, basymptote, ntest, tol)
            [idMax, idSL, idSR, xShockL, xShockR, Iz] = IS_ShockBoundaries(Profile, Soliton, Shock, B13best, X13best, bo, 0, 0, ntest, tol)
            B13best = B13best(idSL:idSR);
            X13best = X13best(idSL:idSR);
            %         Soliton.Type  = 1   ;  % 1 --> 1p soliton, 2 --> 2p soliton, 3 --> KNS soliton
            %         Soliton.lambda = 0.02;  % For a 1p soliton & actual eigenvalue is lambda_actual = lambda*bo
            %         Soliton.sign  = 1  ;  % +1 --> bright 1p soliton, -1 --> dark 1p soliton
            %         Soliton.lr    = 0.5 ;  % Real part of eigenvalue for a 2p soliton --> actual eigenvalue is lambda_actual = (lr + I*li) *bo
            %         Soliton.li    = 0.4 ;  % Imaginary part of eigenvalue for a 1p soliton
            %         Soliton.Locate = -1; % -1 --> place soliton on LHS of shock, +1 --> place soliton on RHS of shock

            %         Shock.Type  =  1    ;  % 1 --> 13 IS, 2 --> 23 IS on lower branch, 3 --> 23 IS on upper branch
            %         Shock.CR    = -0.42 ;  % Compression ratio of 13 IS with -1 < Shock.CR < 0 and CR23 = -1 - IS.CR13
            %         % For Fast Shock 0 < CR < 1
            %         Shock.Iz    = - 1.52; % Target Iz for IS
            %
            %         % theta_L and theta_kbw are only used in the KBW profile
            %         Shock.theta_L = 0*pi; % Normally set at 0 to match a Bfield in the vertical direction at the far LHS (by > 0 and bz = 0)
            %         Shock.theta_kbw = 175/180*pi*(- sign(Shock.Iz) ); % Direction of Bfield downstream of shock (only used in KBW)
            %
            %         % Specify location for shock and ramp
            %         Shock.x     = -0.4395*L;  % location of midpoint of Shock in range: -L --> +L
            %         Shock.xRamp = +0.3469*L;  % location of ramp ... in range: -L --> +L

            % bb --> section of B-field with IntermediatE Shock on LHS preceeded
            % by section with constant B-field = ( bL, 0)
            %
            % x --> array of position values of same length as bb
            %
            % bL ---> as described in bb above
            %
            % basymptote --> it is assumed that the z-component of the B-field
            % goes to this constant asymptotic value to the right of the shock.
            % It is should be evaluated prior to call as basymptote = b_L*sin(theta_kbw)
            %
            %ntest --> number of indices over which the running average of Bz
            %is evaluated


            dX13 = min( X13best(2:end) - X13best(1:end-1) ); % minimum grid spacing needed to match finest resolution of Shock
            Length13 = (X13best(end) - X13best(1)); % length of 1->3 shock
            N13IStmp =floor( 1 + Length13/dX13);      % length of B13
            N13WholeRange = floor(2*L/Length13*N13IStmp);  % # uniform grid points across whole range of 2L given dx = dX13
            N13WholeRange = N13WholeRange + mod(N13WholeRange,2);  % adjusts total # of grid points to be an even number

            XISLeftEdge = xIS - Length13/2;
            XISRightEdge = xIS + Length13/2;

            if XISLeftEdge < (-L)
                XISLeftEdge = -L*.97;
                XISRightEdge = XISLeftEdge + Length13;
                ShockBCmsg1 = sprintf('Given xIS = %4.2f is beyond L = %4.2f, IS shifted to align with left edge of window.',xIS,Length13,XISLeftEdge,-L);
                ShockBCmsg1
            end

            if N13WholeRange <= N % if resolution of shock is more course than default
                N13IS = floor(N/N13WholeRange*N13IStmp); % the larger number of grid points needed for shock in finer resolution
                Nupdate = N;
            else                  % In this case, the shock resolution is finer or equal to default
                                  % thus, the number of grid points needed to resolve the shock remains at N13IStmp
                Nupdate = floor(N13WholeRange);
                N13IS = N13IStmp;
    %             PStr.NCon('N') = Nupdate;                 % PStr not defined, generically, until prior to call to Numerical Method ... overwrite N with updated value
                if mod(Nupdate,2) == 1
                    N = Nupdate+1;
                else
                    N = Nupdate;
                end
            end

            X13uniform = linspace(XISLeftEdge, XISRightEdge,N13IS);
            X13interpolate = linspace(X13best(1), X13best(end),N13IS);

            b13IS = interp1(X13best,B13best,X13interpolate); % This is the 1->3 Shock interpolated to match a uniform grid
            X13uniform = linspace(XISLeftEdge, XISRightEdge,N13IS); % This is the section of the grid where Shock is to be located
            XKBW_IS = linspace(-L, L,Nupdate);              % This is the entire grid from -L to +L with Nupdate => N points
        end
    %     %% Update Prfl Parameters for PrflKBW
    %     if ShockFound
    %         %*********************************************************************
    %         %   Update all profile parameters for KBW profile, a.k.a. PrflKBW, ...
    %         %   before bkbw(PrflKBW,x)
    %         %*********************************************************************    
    %         v = (1 + (Shock.CR)^2 + (Shock.CR) )*alpha*abs(b_L)^2;
    %         CR = -1/2*(1 + sqrt(1 + 4*( v - b_L^2)));
    %         lambda3 = v/(alpha*b_L^2) - 1;
    %         r3 = sqrt(lambda3 + 0.25);
    %         q13 = -1/2 + r3;
    %         q23 = -1/2 - r3;
    %         PrflKBWKeys   = {'title',    'Profile', 'speed', 'CR13', 'CR23', 'Wk'            , 'Wr' , 'theta_L', 'theta_kbw', 'Iz'      , 'theta_R' , 'b_L', 'b_kbw'  , 'b_R','Xkbw'        ,'Xr'   , 'boolShockFit'};
    % %         PrflKBWValues = { Prfltitle,  Profile ,  v     ,  CR13 ,  (-1 - CR13)  ,  L/40           ,  WR  ,  theta_L ,  theta_IS  ,  Iztarget ,  theta_R  ,  b_L ,  b_kbw   ,  b_R , xIS   , xRkbw ,  boolShockFit };
    %         PrflKBWValues = { Prfltitle,  Profile ,  v     ,  q13  ,  q23  ,  L/40/12.2      ,  WR  ,  theta_L ,  theta_IS  ,  Iztarget ,  theta_R  ,  b_L ,  b_kbw   ,  b_R , xIS - 0.3*(xIS + L)  , xRkbw ,  boolShockFit };        
    %         PrflKBW       = containers.Map(PrflKBWKeys,PrflKBWValues);  
    %     end
    %     %figure;plot(XKBW_IS,real(bkbw(PrflKBW,XKBW_IS)));hold 'on';plot(XKBW_IS,imag(bkbw(PrflKBW,XKBW_IS)));hold 'on';plot(XKBW_IS,abs(bkbw(PrflKBW,XKBW_IS)));
        %% Create 1->3 IS Profile on matching background with ramp on RHS for periodic boundary conditions
        if ShockFound
            bKBW_IS = zeros(1,Nupdate);
            idx13IS = 1;
            for i13 = 1:Nupdate
                if XKBW_IS(i13) <  XISLeftEdge
                    bKBW_IS(i13) = b_L;
                elseif (XKBW_IS(i13) >=  XISLeftEdge) && (XKBW_IS(i13) <=  XISRightEdge) &&(idx13IS <= N13IS)
                    bKBW_IS(i13) = b13IS(idx13IS);
                    idx13IS = idx13IS + 1;
                else
                    bKBW_IS(i13) = KBWprofile(PrflKBW, Physics, XKBW_IS(i13), 0);
                end
            end
        end
        %figure;plot(XKBW_IS,real(bKBW_IS));hold 'on';plot(XKBW_IS,imag(bKBW_IS));hold 'on';plot(XKBW_IS,abs(bKBW_IS));
        %% Show Graphs for debugging
        if ShockFound
            if boolDebugProfileGraphs
                %*********************************************************
                %   Graph bkbw as well as B13
                %********************************************************* 

                close all;

                BFieldWindowPosition = [(NMonitors - 1 +0.3) 0.65 0.4 0.3];
                BFieldFigName = 'Magnetic Field Graph';
                BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition);
                BFieldPlot(BFieldFig,b13IS,X13uniform);

                KBWFieldWindowPosition = [(NMonitors - 1 +0.3) 0.35 0.4 0.28];
                KBWFieldFigName = 'KBW Field Graph';
                KBWFieldFig = figure('Name',KBWFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',KBWFieldWindowPosition);
                BFieldPlot(KBWFieldFig,bkbw(PrflKBW,XKBW_IS),XKBW_IS);

                ISShockFieldWindowPosition = [(NMonitors - 1 +0.3) 0.05 0.4 0.28];
                ISShockFieldFigName = 'Hybrid 1->3 Shock over KBW Profile';
                ISShockFieldFig = figure('Name',ISShockFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',ISShockFieldWindowPosition);
                BFieldPlot(ISShockFieldFig,bKBW_IS,XKBW_IS);
                %*********************************************************    
            end
        end    
        %% Update Prfl Parameters for Generated13IS
        if ShockFound
             %*********************************************************************
            %   Update all profile parameters for KBW profile, a.k.a. PrflKBW, ...
            %   before bkbw(PrflKBW,x)
            %*********************************************************************    
            v = (1 + (Shock.CR)^2 + (Shock.CR) )*alpha*abs(b_L)^2;
            CR = -1/2*(1 + sqrt(1 + 4*( v - b_L^2)));
            lambda3 = v/(alpha*b_L^2) - 1;
            r3 = sqrt(lambda3 + 0.25);
            q13 = -1/2 + r3;
            q23 = -1/2 - r3;
            if (ProfileType == 9)
                PrflGenISKeys   = {  'title', 'Profile' , 'speed', 'CR13', 'CR23', 'Wk'            , 'Wr' , 'theta_L', 'theta_kbw', 'Iz_actual', 'Iz_requested', 'theta_R' , 'b_L', 'b_kbw', 'b_R','Xkbw'                ,'Xr'   , 'boolShockFit'};
            else
                PrflGenISKeys   = {  'title', 'Profile2', 'speed', 'CR13', 'CR23', 'Wk'            , 'Wr' , 'theta_L', 'theta_kbw', 'Iz_actual', 'Iz_requested', 'theta_R' , 'b_L', 'b_kbw', 'b_R','Xkbw'                ,'Xr'   , 'boolShockFit'};
            end

                PrflGenISValues = {Prfltitle,  Profile  ,  v          ,  q13  ,  q23  ,  L/40/12.2      ,  WR  ,  theta_L ,  theta_IS  ,  Iztarget  , Shock.Iz      ,  theta_R  ,  b_L ,  b_kbw ,  b_R , xIS - 0.3*(xIS + L)  , xRkbw ,  boolShockFit };        
            PrflGeneratedIS       = containers.Map(PrflGenISKeys,PrflGenISValues);   
            %*********************************************************************  
            %   Update NumericsGeneratedIS ... if N has changed, it will need
            %   to be changed everywhere.
            %*********************************************************************
    %         NumericsKeys = {'title', 'N', 'L', 'dx', 'Courant', 'dt', 'timestepsPerFrame', 'TimestepsAttempted', 'TotalTime'};
    %         NumericsValues = {Ntitle, N, L, dx, Courant, dt, timestepsPerFrame, TimestepsAttempted, TotalTime};
    %         Numerics = containers.Map(NumericsKeys, NumericsValues);
            Numerics('N') = N;

            %*********************************************************************  
            %   Update GraphicsGeneratedIS ...
            %*********************************************************************
            zmax = 1.5*max(abs(bKBW_IS));
            zmin = -zmax;
            GraphicsValues = {Gtitle, exact, zmax,zmin, NMonitors};
            GraphicsGeneratedIS = containers.Map(GraphicsKeys, GraphicsValues);    

            %*********************************************************************
            % Primary Result of GeneratedIS:
            %  PrflGeneratedIS
            %  GenISprofile: The generated Intermediate Shock consistent with parameterw in PrflGeneratedIS
            %*********************************************************************
            GenISprofile = @(Prfl, Physics, x,t) bKBW_IS;

            % If this is called as a stand alone, then define Prfl and profile
            % ...
            if (ProfileType == 9)
                Prfl     = PrflGeneratedIS;
                Graphics = GraphicsGeneratedIS;
                PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
                    'PrflCon', Prfl, 'GCon',Graphics,'FCon',Folder, 'DrvrCon',Driver);
                profile = @(Prfl, Physics, x,t) bKBW_IS;
            end
            Proceed = 1;
        else
            msgBG
            Proceed = 0;
        end
        % Results in bKBW_IS with parameters defined in PrflKBW
        % Iztarget becomes the shock Iz actually found
     end
    %% ProfileType #10     Continuation of previous trial
    if (ProfileType == 10)
        %% Trial, Run and Line Number to Continue
        TrialNumberOld = '60';     % Trial # to be continued
        RunNumberOld = '4';       % Run # to be continued
        LineNumberOld = 100;       % Line # of TrialNumberOld, RunNumberOld to use as the initial profile
        %If LineNumberOld > (the # of lines%recorded), then LineNumberOld != last line # recorded

        % The new TrialNumber, RunNumber_start, RunNumber_end must be declared
        % above in "Trial to Create or Overwrite"
        RunNumber_end   = RunNumber; % Sorry, multiple trial features have not been implemented.
        RunNumberNew    = num2str(RunNumber);
        %% Read PSTrOld and Bfield for LineNumberOld    
        % read last archived line if LineNumberOld is not present
        Verbal = 0;
        [PStrOld, B4NewStart, LineNumberOldRead, MSG1, flag1] = Continuation_Read_Bfield_Line(ArchiveBaseOld, ActiveArchiveBase, TrialNumberOld, RunNumberOld, LineNumberOld, TrialNumber, RunNumberNew, Verbal);
        %% Update and Save various parameters to NCon, TCon and PrflCon
        PStr = PStrOld;    
        Prfl = PStr.PrflCon;
        NCon = PStr.NCon;
        TCon = PStr.TCon;

        % Used to interpolate Bfield
        N1 = NCon('N');
        L1 = NCon('L');

        % Save Trial, Run and LineNumber to be continued to Trial Identifiers
        TrialIDOld   = strcat(TrialNumberOld,'_',RunNumberOld,'_Line_',int2str(LineNumberOldRead));
        TrialID   = strcat(TrialNumber,'_',int2str(RunNumber)); % "#-#' format
        TCon('OriginalTrial') = TrialIDOld;
        TCon('TrialID')       = TrialID;
        TCon('TrialNumber')   = TrialNumber;
        TCon('RunNumber')     = RunNumberNew;    

        % Update TotalTime, TimestepStart, LineNumberStart and save to Numerics Parameters
        dtOld              = PStrOld.NCon('dt');
        TimestepsAttempted = round(TotalTime/dtOld) + 1; % # timesteps to achieve in continuation
        timestepsPerFrame  = PStrOld.NCon('timestepsPerFrame');
        TrialKeys = TCon.keys; 
        TrialValues = TCon.values;
        TCon = containers.Map(TrialKeys, TrialValues);

        NCon('TotalTime')          = TotalTime;
        NCon('TimestepStart')      = timestepsPerFrame*(LineNumberOldRead - 1);
        NCon('LineNumberStart')    = LineNumberOldRead;
        NCon('TimestepsAttempted') = TimestepsAttempted;
        NumericsKeys   = NCon.keys;
        NumericsValues = NCon.values; 

        %% Update Prfl with 'Continuation'
        PrflKeys   = Prfl.keys;
        PrflValues = Prfl.values;    
        ContinueAgain = strfind(PrflKeys,'ProfileNew');
        CopyofCopy = any( vertcat( ContinueAgain{:} ) );
        if CopyofCopy
            CopyMsg = ['I am just a copy of a copy of a copy',newline,...
                'Everything I say has come before',newline,...
                'Assembled into something, into something, into something',newline,...
                'I don''t know for certain anymore',newline,...
                'https://genius.com/Nine-inch-nails-copy-of-a-lyrics'];    
            CopyMsg
        end
        Prfl('ProfileNew') = ['Continuation of ',Prfl('Profile')];
        PrflKeys   = Prfl.keys;
        PrflValues = Prfl.values;

        %% Interpolate Bfield to match original grid
        X4uniform    = linspace(-L1, +L1, length(B4NewStart));
        Xinterpolate = linspace(-L1, +L1, N1);   
        BNewStart = interp1(X4uniform, B4NewStart, Xinterpolate); % B-field from old trial to be passed to DNLS_v4 as start of extended trial

        %% Update NCon, PrflCon, TCon, Graphics, Driver and store in PStrContinuation
        NCon     = containers.Map(NumericsKeys, NumericsValues);
        PrflCon  = containers.Map(PrflKeys,PrflValues);
        % TCon updated directly above
        zmax = 1.5*max(abs(BNewStart));
        zmin = -zmax;
        GraphicsValues = {Gtitle, exact, zmax,zmin, NMonitors};
        Graphics = containers.Map(GraphicsKeys, GraphicsValues);

        DriverKeys   = {'title' ,     'IsDriven'};
        DriverValues = { Drivertitle , IsDriven };
        Driver = containers.Map(DriverKeys, DriverValues);

    %     FolderValues = {Ftitle,  ArchiveBase ,  TrialFolderAbs, BfieldFolderAbs, B_field_Tr_Rn_txt, B_field_Tr_mat};
    %     FolderKeys = {'title' , 'ArchiveBase', 'TrialFolder'  , 'BfieldFolder' ,   'Bfieldfile'   , 'BfieldMatfile'  };
    %     Folder = containers.Map(FolderKeys,FolderValues);
    %     FCon

        PStrContinuation = struct('TCon',TCon,'NCon',NCon,'PCon',Physics,...
            'PrflCon', PrflCon, 'GCon',Graphics, 'DrvrCon',Driver);

        Proceed = flag1;


        % ToDo:
        %       
        %       X1) Determine variables to update ... and update them.
        %       X2) Determine variables that are necessary to pass to DNLS_v4
        %       X3) Create extended and modified PStr and archive
        %       4) Read all Bfield lines from old trial
        %       5) Archive (LineNumber - 1) of these to txt and mat files
        %           --> Logistically, it is more straightforward to write these
        %           lines to the Bfile from DNLS_v4 ... just before
        %           timestepping starts.
        %       X6) Prepare Bfield LineNumberN interpolated to N to pass to
        %       DNLS_V4 through 'Profile'
        %       7) Bfield lines up to lineN -1 to txt and matfile
        %       X8) While PStr is updated here, prevent PStr from being
        %       overwritten during Archiving features prior to Numerical Method
        %       X9) Modify DNLS_v4 to check for 'Continuation' and to use values
        %       of timesteps, linenum and TimestepAttempted

    end
    %% ProfileType #11     Hybrid: Generated Intermediate Shock vs Soliton
    if (ProfileType == 11) &&  ProceedKBW == 1
    Profile = 'Hybrid Shock-Soliton';

        %% Prepare Soliton Selection
        %***********************************************************    
        % Prepare Soliton: 
        %  --> Set up Prfl2, the container with ...
        %      finalized soliton type and sign (if 1p) & 'eigenfactors':Hybrid.Locate
        %      eigenvalue = (eigenfactor)*background field strength
        %  --> determine Lsoliton from graph
        %  --> Other Prfl2 parameters, such b_soliton, will be set
        %      once the background field strength is determined based on
        %      b_L, the field strength to the left of the shock, and
        %      the shock's compression ratio and, eventually, the field
        %      rotation across the shock. Currently, it is assumed to by 180 degrees
        %***********************************************************     
        %     Hybrid.Type   =  1; %  1 --> Archived Trial vs. soliton, 2 --> Generated IS vs. soliton
        %     Hybrid.Locate = -1; % -1 --> place soliton on LHS of shock, +1 --> place soliton on RHS of shock
        % b_soliton temporary until splice location is set
    %     Shock.Type  =  1    ;  % 1 --> 13 IS, 2 --> 23 IS on lower branch, 3 --> 23 IS on upper branch
    %     Shock.CR    = -0.42 ;  % Compression ratio of 13 IS with -0.5 < IS.CR13 < 0 and CR23 = -1 - IS.CR13
    %                            % For Fast Shock 0 < CR < 1
    %     Shock.Iz    = - 1.3; % Target Iz for IS
    %     
    %     % theta_L and theta_kbw are only used in the KBW profile
    %     Shock.theta_L = 0*pi; % Normally set at 0 to match a Bfield in the vertical direction at the far LHS (by > 0 and bz = 0)
    %     Shock.theta_kbw = 175/180*pi*(- sign(Shock.Iz) ); % Direction of Bfield downstream of shock (only used in KBW)
    %     
    %     % Specify location for shock and ramp
    %     Shock.x     = -0.4395*L;  % location of midpoint of Shock in range: -L --> +L
    %     Shock.xRamp = +0.3469*L;  % location of ramp ... in range: -L --> +L 
    if (Soliton.Locate == -1)
        b_soliton = bo;
    elseif Soliton.Locate == 1
        b_soliton = abs(Shock.CR)*bo;
    end

    %     Soliton.Type  = 1   ;  % 1 --> 1p soliton, 2 --> 2p soliton, 3 --> KNS soliton
    %     Soliton.lambda = 0.02;  % For a 1p soliton & actual eigenvalue is lambda_actual = lambda*bo
    %     Soliton.sign  = 1  ;  % +1 --> bright 1p soliton, -1 --> dark 1p soliton
    %     Soliton.lr    = 0.5 ;  % Real part of eigenvalue for a 2p soliton --> actual eigenvalue is lambda_actual = (lr + I*li) *bo
    %     Soliton.li    = 0.4 ;  % Imaginary part of eigenvalue for a 1p soliton
    %% Set UP PrflSoliton and profileSoliton(PrflSoliton, Physics, x, t)
    if Soliton.Type == 1
        % 1p-Soliton parameters
        lambda   = Soliton.lambda; % For 1p_soliton: lambda = LambdaFactor*(local field strength) with % 0 < LambdaFactor < 1
        sign1 = Soliton.sign;       % sign = +1 ... 'bright soliton'; sign = -1 ... for a 'dark soliton'
        Profile1 = '1psoliton';

        %*********************************************************
        %   update PrflSoliton (the soliton parameters) and
        %   profileSoliton (the B-Field for the soliton)
        %*********************************************************
        vsol = alpha*b_soliton^2*(1 + 2*lambda^2);  % the soliton speed ... to be updated with change
        % to lambda or bl & just prior to update of Prfl

        PrflOnePSolitonKeys   = {'title'   , 'Profile1', 'Solitonspeed', 'lambda',  'sign'    , 'b1p'       , 'Lscale'                 };
        PrflOnePSolitonValues = { Prfltitle,  Profile1 ,  vsol         ,  lambda ,   sign1 ,  b_soliton  ,  1/(lambda*b_soliton^2)  };
        PrflSoliton           = containers.Map(PrflOnePSolitonKeys,PrflOnePSolitonValues);

        profileSoliton = @(PrflSoliton, Physics, x, t) (onePsoliton(PrflSoliton, Physics, x, t));

        debug1pgraph = 1;
        if debug1pgraph
            close all;
            X1p = linspace(-L,L,N);
            BFieldWindowPosition = [(NMonitors - 1 +0.3) 0.75 0.4 0.2];
            BFieldFigName = sprintf('%s Graph', 'onePsoliton');
            BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition,'visible','off');
            BFieldPlot(BFieldFig,onePsoliton(Prfl1p,Physics,X1p,0),X1p);
        end

    %     debug1pgraph = 1;
        if debug1pgraph
            X1p = linspace(-L,L,N);
            BFieldWindowPosition = [(NMonitors - 1 +0.3) 0.5 0.4 0.2];
            BFieldFigName = sprintf('%s Graph', 'profileSoliton with Prfl1p');
            BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition,'visible','off');
            BFieldPlot(BFieldFig,profileSoliton(Prfl1p,Physics,X1p,0),X1p);
        end

    %     debug1pgraph = 1;
        if debug1pgraph
            X1p = linspace(-L,L,N);
            BFieldWindowPosition = [(NMonitors - 1 +0.3) 0.15 0.4 0.3];
            BFieldFigName = sprintf('%s Graph', 'profileSoliton with PrflSoliton');
            BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition,'visible','off');
            BFieldPlot(BFieldFig,profileSoliton(PrflSoliton,Physics,X1p,0),X1p);
        end


        %*********************************************************
        prfmsg = sprintf('Finalized 1psoliton profile with %s = (%4.4f)bo with bo = %4.4f and sign = %4.0f',lambdaChar,lambda,b_soliton, sign1);
        prfmsg
    elseif Soliton.Type == 2
        % 2p-Soliton parameters
        lr = Soliton.lr;   % set here s.t. actual Re(lambda) = lr*(local field strength)
        li = Soliton.li;   % set here s.t. actual Im(lambda) = li*(local field strength)
        Profile1 = '2psoliton';

        %*********************************************************
        %   update PrflSoliton (the soliton parameters) and
        %   profileSoliton (the B-Field for the soliton)
        %*********************************************************
        vsol = alpha*b_soliton^2*(1 + 2*(lr^2 - li^2));  % the soliton speed ... to be updated with change
        % to lambda or bl & just prior to update of Prfl

        PrflTwoPSolitonKeys   = {'title'   , 'Profile1', 'speed'   , 'lr','li' , 'b2p'      , 'Lscale'};
        PrflTwoPSolitonValues = { Prfltitle,  Profile1 ,  vsol     ,  lr , li  ,  b_soliton ,  0.2*L  };
        PrflSoliton           = containers.Map(PrflTwoPSolitonKeys,PrflTwoPSolitonValues);

        PrflTwoPSolitonKeys   = {'title'   , 'Profile1' , 'Solitonspeed', 'lr','li' , 'b2p'      , 'Lscale'                      };
        PrflTwoPSolitonValues = { Prfltitle,  Profile1  ,  vsol         ,  lr , li  ,  b_soliton ,  Lscale(PrflSoliton, Physics) };
        PrflSoliton           = containers.Map(PrflTwoPSolitonKeys,PrflTwoPSolitonValues);

        profileSoliton = @(PrflSoliton, Physics,x,t) twopSoliton(PrflSoliton, Physics,x,t);
        %*********************************************************
        prfmsg = sprintf('Finalized 2psoliton profile with %s = ( %4.4f + i(%4.4f) )bo with bo = %4.4f',lambdaChar,lr, li, b_soliton);
        prfmsg
    end
    % Result is: PrflSoliton, whose values can be adjusted below as needed, to create bsol ... the soliton B-field
    %% Establish Soliton Boundaries
    % Profile specific values
    NTestSoliton = 10; Tol_Soliton = 0.001*bo; b_ref = 1; % Potentially useful in later development to allow user to adjust tolerance in selecting boundaries
    % to adjust tolerance within interactive while loop
    Prfl_Length_key = 'Lsoliton';
    Xprofile = linspace(-L,L,N);                     % Assumption that profile is defined on entire range
    bsol = profileSoliton(PrflSoliton, Physics, Xprofile, 0);  % Define the soliton B-field with boundaries to identify
    bfit = bsol;                                       % Default profile name to pass to boundary approval routine

    %Define figure for B-field soliton graph used to verify soliton boundaries
    BFieldWindowPosition = [(NMonitors - 1 +0.3) 0.65 0.4 0.3];
    BFieldFigName = sprintf('%s Graph', PrflSoliton('Profile1'));
    B1flag = 1;                                       % Changed to another value in while loop once user approves boundaries
    idL1 = 1;idR1 = N;                                % Assumption that profile is defined on entire range
    xL = -L; xR = +L;                                  % Assumption that profile is defined on entire range
    zmax = 1.5*max( abs(bfit) );                      % Determine y-scale for graphing
    zmin = -zmax;
    comment ='Left click to select region.';
    FirstPass = 1;

        while B1flag == 1
            BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition,'visible','off');
            clf;
    %         axis([Xprofile(idL1) Xprofile(idR1) zmin zmax]);
            TmpProfile = [num2str(Soliton.Type),'p-Soliton'];
            [idRelMax, idxShockL, idxShockR, XShockL, XShockR, Izcalc] = IS_ShockBoundaries(TmpProfile, Soliton, Shock, bfit, Xprofile, b_soliton, 0, 0, NTestSoliton, Tol_Soliton);
            DeltaIdShock = idxShockR - idxShockL;
            idL1 = idxShockL - floor(DeltaIdShock);
            if idL1 < 1
                idL1 = 1;
            end

            idR1 = idxShockR + floor(DeltaIdShock);
            if idR1 > N
                idR1 = N;
            end


            axis([Xprofile(idL1) Xprofile(idR1) zmin zmax]);
            IS_Shock_Profile_Label(BFieldFig, TmpProfile, 0, zmax,XShockL, XShockR, Xprofile(idL1:idR1), bfit(idL1:idR1), comment);
            BFieldPlot(BFieldFig,bfit(idL1:idR1),Xprofile(idL1:idR1));
            [B1flag, FirstPass, Prfl_fit, xL_profile, xR_profile, L_profile, idL1, idR1] = Boundary_Approval(BFieldFig, TmpProfile, PrflSoliton, Prfl_Length_key, FirstPass, Xprofile, XShockL, XShockR);
            close(BFieldFig);
        end
        xL_soliton = xL_profile;
        xR_soliton = xR_profile;
        Lsoliton   = L_profile;
        idLsol = idL1;
        idRsol = idR1;
        % Result is: bsol, xL_soliton, xR_soliton, idLsol, idRsol and Lsoliton
        % along with Prfl2 (with Lsoliton as a key/value pair)

        %% Prepare Shock Selection

        %     {bKBW_IS,XKBW_IS, Iztarget} returned from Generate13IS

        %*********************************************************************

        Profile2 = sprintf('Shock');
        if Shock.Type == 1
            Profile2 = sprintf('1->3 Shock');
            vIS = (1 + CR^2 + CR )*alpha*abs(PrflKBW('b_L'))^2;
        elseif Shock.Type > 1
            Profile2 = sprintf('2->3 Shock');
            vIS = (1 + CR^2 + CR )*alpha*abs(PrflKBW('b_L'))^2;
        end

    %     PrflKBWKeys   = {'title',    'Profile', 'speed'  , 'CR13', 'CR23'        , 'Wk'                 , 'Wr' , 'theta_L', 'theta_kbw', 'Iz'      , 'theta_R' , 'b_L', 'b_kbw' , 'b_R','Xkbw'          ,'Xr'   , 'boolShockFit', 'SolitonLocation'};
    % %     PrflKBWValues = { Prfltitle,  Profile ,  vIS     ,  CR13 ,  (-1 - CR13)  ,  1/(2*12.2)*Length13 ,  WR  ,  theta_L ,  theta_IS  ,  Iztarget ,  theta_R  ,  b_L ,  b_kbw  ,  b_R , xIS         , xRkbw ,  boolShockFit ,   HybridIS.Locate};
    %     PrflKBWValues = { Prfltitle,  Profile ,  vIS     ,  CR13 ,  (-1 - CR13)  ,  1/(2*12.2)*Length13 ,  WR  ,  theta_L ,  theta_IS  ,  Iztarget ,  theta_R  ,  b_L ,  b_kbw  ,  b_R , xIS - Length13 , xRkbw ,  boolShockFit ,   HybridIS.Locate};
    %     Prfl1       = containers.Map(PrflKBWKeys,PrflKBWValues);   
        %********************************************************************  

        %% Establish Shock Boundaries
        % Profile specific values
        NTestShock = 10; Tol_Shock = 0.001; b_ref = 1; % Potentially useful in later development to allow user to adjust tolerance in selecting boundaries
                                                       % to adjust tolerance within interactive while loop
        Prfl_Length_key = 'LShock';
        Prfl_fit = PrflGeneratedIS;         % Prfl1 is the Prfl container for the Shock
        bfit = GenISprofile(Prfl_fit, Physics, Xprofile, 0);           % bfit is default profile name to pass to boundary approval routine


        %Define figure for B-field shock graph used to verify shock boundaries
        BFieldWindowPosition = [(NMonitors - 1 +0.3) 0.65 0.4 0.3]; 
        BFieldFigName = sprintf('%s Graph', Prfl_fit('Profile2'));  
        B1flag = 1;                                       % Changed to another value in while loop once user approves boundaries
        idL1 = 1;idR1 = N;                                % Assumption that profile is defined on entire range
        xL = -L;xR = +L;                                  % Assumption that profile is defined on entire range
        Xprofile = linspace(xL,xR,N);                     % Assumption that profile is defined on entire range 
        zmax = 1.5*max( abs(bfit) );                      % Determine y-scale for graphing
        zmin = -zmax;
        comment ='Left click to select region.';
        FirstPass = 1;

        while B1flag == 1
            BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition,'visible','off');
            clf;
    %         axis([Xprofile(idL1) Xprofile(idR1) zmin zmax]);
    %         [idRelMax, idxShockL, idxShockR, XShockL, XShockR, Izcalc] = IS_ShockBoundaries(Profile, bfit, Xprofile, b_ref, 0, 0, NTestShock, Tol_Shock);
            XShockL = XISLeftEdge;
            XShockR = XISRightEdge;
    %         DeltaIdShock = idxShockR - idxShockL;
    %         idL1 = idxShockL - floor(DeltaIdShock/2);
    %         if idL1 < 1
    %             idL1 = 1;
    %         end
    %        
    %         idR1 = idxShockR + floor(DeltaIdShock/2);
    %         if idR1 > N
    %             idR1 = N;
    %         end

            axis([Xprofile(idL1) Xprofile(idR1) zmin zmax]);
            IS_Shock_Profile_Label(BFieldFig, Profile,0, zmax, XShockL, XShockR, Xprofile(idL1:idR1), bfit(idL1:idR1), comment);
            BFieldPlot(BFieldFig,bfit(idL1:idR1),Xprofile(idL1:idR1));
            [B1flag, FirstPass, Prfl_fit, xL_profile, xR_profile, L_profile, idL1, idR1] = Boundary_Approval(BFieldFig, Profile, Prfl_fit, 'Lsoliton', FirstPass, Xprofile, XShockL, XShockR);
            close(BFieldFig);
        end
        XShockL = xL_profile;
        XShockR = xR_profile;
        LShock   = L_profile;
        tmpL1 = abs(Xprofile - xR_profile);
        [val, idxShockR] = min(tmpL1);
        tmpL1 = abs(Xprofile - xL_profile);
        [val, idxShockL] = min(tmpL1);
        % Result is: bsol, xL_soliton, xR_soliton, idLsol, idRsol and Lsoliton
        % along with Prfl2 (with Lsoliton as a key/value pair)   

        %% Combine Shock and Soliton
        axis([-L, L, zmin, zmax]);
        Zmax = 1.5*max(max(abs(bKBW_IS)),max(abs(bsol)));
        Zmin = -Zmax;
        BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition,'visible','off');
            clf;
        if Soliton.Locate == -1

            idxRsoliton = idxShockL;
            xLsoliton = XShockL - Lsoliton;
            xRsoliton = XShockL;
            tmpL1 = abs(Xprofile - xLsoliton);
            [val, idxLsoliton] = min(tmpL1);

        elseif Soliton.Locate == +1

            idxLsoliton = idxShockR;
            xRsoliton = XShockR + Lsoliton;
            xLsoliton = XShockR;
            tmpL1 = abs(Xprofile - xRsoliton);
            [val, idxRsoliton] = min(tmpL1);

        end

        % Suggestions for optimal locations of Shock and Ramp
        % Given optimal locations: 
        %  --> estimate of minimum required run time for Shock-Soliton collision
        %  --> estimate of maximum time of trial validity before the shock or
        %      the soliton will collide with the ramp

            %% Ramp width

            % the ramp transition on the RHS goes as tanh(x/WR)
            % for ramp on RHS to be with epsilong of RHS boundary conditions,then ...
            % tanh(Lramp/2/WR) = 1 - epsilon, with x = 0 set to match the ramp position, 
            % 'Lramp/2' is used as the ramp extends a distance Lramp/2 to both the left and the right
            % As tanh(x/WR) = (exp(+x/WR) - exp(-x/WR)/(exp(+x/WR) + exp(-x/WR), then with y = exp(+x/WR) ...
            % (y - 1/y)/(y + 1/y) = 1 - e --> y = ( (2 - e)/e )^(1/2) 
            % --> exp(Lramp/2/WR) = ( (2 - e)/e )^(1/2) 
            % --> Lramp/2/WR = Ln( ( (2 - e)/e )^(1/2) ) =1/2 Ln( ( (2 - e)/e  ) = 1/2 Ln(  2/e ) to first order in e << 1
            % --> Lramp = WR*Ln(2/e) = 12.2*WR for e = 10^(-5)
            Lramp = 12.2*WR;

            % ProfileError = 0 ... all is good
            % ProfileError = 1 ... shock + soliton + ramp are too wide to fit in region (-L,L)
            % ProfileError = 2 ... shock-soliton extends beyond LHS of (-L, L)
            % ProfileError = 3 ... ramp extends beyond RHS of (-L, L)
            % ProfileError = 4 ... shock-soliton overlaps LHS of ramp
            % ProfileError = 5 ... shock-soliton might not collide
            PrfErrMsg = "Everything is going just great ... so far.";
            ProfileError = 0;


            LPrfl = Lsoliton + LShock + Lramp; % width of shock, soliton and ramp ... if not overlapping each other

            %Left and Right edges of Ramp
            XrampL = Shock.xRamp - Lramp/2;
            XrampR = Shock.xRamp + Lramp/2;

            % distance from left edge of ramp to right edge of 'Shock-Soliton' profile
            DXRight = XrampL - max(xRsoliton, XShockR); % if DXRight < 0 --> LHS of ramp and RHS of profile overlap ... BAD!

            % distance from right edge of ramp to RHS = L - XrampL
            % distance from left edge of 'shock-soliton' profile = min(XLsoliton, XShockL) - ( -L)
            DXLeft = ( L - XrampR ) + ( min(xLsoliton, XShockL ) - (-L) ); % by construction this should always be positive

            DXfree = 2*L - LPrfl; 
            % Total 'flat line' distance. 
            % If this is split in half between RHS of profile to ramp and LHS of profile and ramp, 
            % then the ramp's rarefactive wave can travel approximately
            % DXfree/2 before colliding either from the right, or the left,
            % with the profile, thus, destroying the 'physical' boundary
            % conditions required for the simulation.
            % The rarefactive wave appears to travel at a speed approximated
            % by: VRW = alpha*(bL^2 - bshock^2), where bL is the magnetic field
            % boundary conditions to the far left and bshock is the magnitude
            % of the magnetic field just to the right of the shock.
    %         VRW = alpha*(b_L^2 - b_kbw^2);
            VRWFoot = 0; % naive, but vaguely sane assumption for speed of foot of ramp
            VRW = alpha*(bo^2 + b_kbw^2 + bo*b_kbw); % for bL = 1 and bkbw = 0.25, VRW observed to be about 3x relative speed between Shock and foot of ramp

            PrfWidth = LPrfl;
            CollisionDistance = (2*L - PrfWidth  )/2; 
            % Returns distance between ramp to the right, or to the left of the 
            % shock-soliton ... assuming the profile is spaced right in the
            % middle.
    %         EstFailTime = CollisionDistance/(0.4);
    %         Ladvised = PrfWidth/2 + .2*400;

            impactSpeed = vIS - vsol; % assuming impact speed remains constant prior to collision and throughout collision, because what else can be done?
             if Soliton.Locate == -1            
                impactSpeed = vsol - vIS;          
            elseif Soliton.Locate == +1            
                impactSpeed = vIS - vsol; 
            end
            impactDistance = LShock + Lsoliton; % assuming collision completion occurs when the furthest edge of the soliton has had time to move past the shock. Assuming oh so many things.
            estImpactTime = impactDistance/impactSpeed;

    %         xRampOptimal = L - CollisionDistance;
    %         xISoptimal = -L + CollisionDistance;

    %         Lnew = 1.1*1/2*( LPrfl + (vIS + VRW)*estImpactTime );
            % Insert Dialog box here to determine if user
            % --> wants to runtime < collision time?
            % --> wants runtime > collision time?
            % Based on response, set TBCFail = 1.1*runtime ... with runtime
            % established by boss
             % Construct a questdlg with three options
             StubbornTime = 0;
             if TotalTime < 1.1*estImpactTime
                 choice = questdlg('Runtime is less than estimated duration of collision:', ...
                     'Keep runtime as is, or increase to observe collision?', ...
                     'Keep','Increase', 'No');
                 % Handle response
                 switch choice
                     case 'Keep'
                         disp([choice 'Keeping runtime as is.'])
                         StubbornTime = 1;
                     case 'Increase'
                         disp([choice 'Increasing runtime enough to see collision.'])
                         TotalTime = 1.1*estImpactTime;
                     case 'No'
                         disp([choice ' I''m out of here!'])
                 end
             else
                 choice = questdlg('Runtime is more than estimated duration of collision:', ...
                     'Keep runtime as is, or decrease to minimum needed to observe collision?', ...
                     'Keep','Decrease', 'No');
                 % Handle response
                 switch choice
                     case 'Keep'
                         disp([choice 'Keeping runtime as is.'])
                     case 'Decrease'
                         disp([choice 'Decreasing runtime to just enough to see collision.'])
                         TotalTime = 1.1*estImpactTime;
                     case 'No'
                         disp([choice ' I''m out of here!'])
                 end

             end


            TBCFail = 1.1*estImpactTime; % set so that boundary conditions survive beyond duration of collision
            runtime = max(TotalTime, TBCFail); % just in case boss wants simulation to run well beyond duration of collision

            % Lnew, below, will allow boundary conditions to remain valid through 
            % collision and beyond that if required by TotalTime
            Lnew = 1/2*(LPrfl + (vIS - alpha*b_kbw^2 + VRW)*runtime); 

            if (0.8*Lnew > L) && ProfileError == 0 % boundary conditions will fail prior to collision or TotalTime with given value of L

                PrfErrMsg = "Boundary conditions will fail prior to either collision or TotalTime with given value of L."
                choice = questdlg('L is too small and BC''s may fail too soon:', ...
                     'Keep doomed L as is or restart and Adjust values?',...
                     'Keep','Adjust', 'No');
                 % Handle response
                 switch choice
                     case 'Keep'
                         disp([choice 'Keeping L as is.'])
                     case 'Adjust'
                         disp([choice 'Advice follows ...'])
                         ProfileError = 1;
                     case 'No'
                         disp([choice ' I''m out of here!'])
                         ProfileError = 1;
                 end

                % The total length of the simulation should be at least be:
                %   Lsoliton + LShock + Lramp
                % However, a rarefactive wave will travel to the right and
                % eventually collide with the LHS of the shock-soliton thus
                % destroying the physical boundary conditions. The distance
                % this wave will travel is:
                %                 DXleft = VRW*runtime
                % Additionally, the shock (assumed to be the fastest wave in
                % the simulation) will need to be able to travel a distance
                %                 DXright = vIS*estImpactTime
                % The minimum simulation distance will be the sum of these
                % distances:
                %    2*L = LPrfl + (vIS + VRW)*runtime
                % There is a 'safety factor' of 1.2 included in the revised,
                % Lnew defined above (or 1.1 for value of L).


                msgL2small = sprintf('Suggested value for L = 1/2*( Lshock + Lsoliton + Lramp + (vIS + VRW)*runtime ) = %4.4f ',Lnew);
                msgL2small

                [XrampNew, XISNew, TBCFail, Tcollision] = HybridOptimize(Lnew, LPrfl, Lramp, Lsoliton, LShock, Soliton, vsol, vIS, VRW, VRWFoot);

                if TBCFail < Tcollision
                    Lnew = 1/2*(LPrfl + (vIS - alpha*b_kbw^2 + VRW)*runtime); 
                    [XrampNew, XISNew, TBCFail, Tcollision] = HybridOptimize(Lnew, LPrfl, Lramp, Lsoliton, LShock, Soliton, vsol, vIS, VRW, VRWFoot);
                end

                Tsuggest = Tcollision;
                msgL2smallTime = sprintf('Suggested TotalTime = 1.2*(Lsoliton + LShock)/impactspeed = %4.4f ',Tsuggest);
                msgL2smallTime           

                msgLrampNew = sprintf('Suggest that ramp be placed at Xramp = (%4.4f)Lnew ',XrampNew/Lnew);
                msgLrampNew

                msgISNew = sprintf('Suggest that shock be placed at XIS = (%4.4f)Lnew ',XISNew/Lnew);
                msgISNew
            end

             % Check that shock-soliton aren't too far to left
            if (min(xLsoliton, XShockL ) < (-L)) && ProfileError == 0
                ProfileError = 2;
                PrfErrMsg = "LHS of profile extends beyond known universe (shift shock-soliton to the right, or increase L).";
            end

            % Check that ramp, for whatever reason, is not too far to the right, like how could that even happen? :-O
            if (XrampR > L) && ProfileError == 0
                ProfileError = 3;
                PrfErrMsg = "RHS of Ramp extends beyond known universe (shift ramp to left, or increase L).";
            end    

            % Check that shock-soliton does not overlap ramp
            if (max(xRsoliton, XShockR) > XrampL) && ProfileError == 0
                ProfileError = 4;
                PrfErrMsg = "Profile overlaps LHS of Ramp (shift ramp to right, or shift shock-soliton to the left)."
                Lnew = L;
                [XrampNew, XISNew, TBCFail, Tcollision] = HybridOptimize(Lnew, LPrfl, Lramp, Lsoliton, LShock, Soliton, vsol, vIS, VRW, VRWFoot);

                msgLrampNew = sprintf('Suggest that ramp be placed at Xramp = (%4.4f)L ',XrampNew/L);
                msgLrampNew

                msgISNew = sprintf('Suggest that shock be placed at XIS = (%4.4f)L ',XISnew/L);
                msgISNew

            end

            % Compare estimated time for completion of collision vs. est. time
            % for boundary condition failure
            TcollisionEnd = 1.2*(LShock + Lsoliton)/impactSpeed;
            TBCFailure = DXfree/2/VRW;

            if (estImpactTime > TBCFail) && ProfileError == 0 && StubbornTime ==0
                ProfileError = 6;
                PrfErrMsg = "Boundary conditions will fail before collision is completed."
                Lnew = 1/2*(LPrfl + (vIS + VRW)*runtime);
                msgNotEnoughSpace = sprintf('Need to increase L to Lnew = L + (Lsoliton + LShock)  = %4.4f ',Lnew);
                msgNotEnoughSpace

                msgTime = sprintf('Suggest setting TotalTime = 1.2*(Lsoliton + LShock)/impactSpeed  = %4.4f ',1.2*(Lsoliton + LShock)/impactSpeed);
                msgTime

                [XrampNew, XISNew, TBCFail, Tcollision] = HybridOptimize(Lnew, LPrfl, Lramp, Lsoliton, LShock, Soliton, vsol, vIS, VRW, VRWFoot);

                msgLrampNew = sprintf('Suggest that ramp be placed at Xramp = (%4.4f)Lnew ',XrampNew/Lnew);
                msgLrampNew

                msgISNew = sprintf('Suggest that shock be placed at XIS = (%4.4f)Lnew ',XISNew/Lnew);
                msgISNew
            end




            % Check if shock-soliton initial speeds are consistent with a
            % collision. The soliton speed may change with damping to
            % eventually allow a collision, so this is not an intrinsically
            % fatal problem.
            if Soliton.Locate == -1            
                if impactSpeed < 0 % impactSpeed = vIS - vsol;
                    ProfileError = 0;
                    PrfErrMsg = "Collision uncertain: Soliton is to the left of the, initially, faster shock and might not catch it.";
                end                
            elseif Soliton.Locate == +1            
                if impactSpeed < 0 % impactSpeed = vIS - vsol;
                    ProfileError = 0;
                    PrfErrMsg = "Collision uncertain: Soliton is to the right of the, initially, slower shock.";
                end    
            end




    % 
    %     msgIS = sprintf('Optimal shock location at %4.4f or (%4.4f)L',xISoptimal,xISoptimal/L);
    %     msgIS
    %     msgRmp = sprintf('Optimal ramp location at %4.4f or (%4.4f)L',xRampOptimal,xRampOptimal/L);
    %     msgRmp
    %     msgTimetoImpact = sprintf('Est. time one Shock-Solition collision %4.4f',estImpactTime);
    %     msgTimetoImpact
    %     msgTimetoFailure = sprintf('Est. time of trial validity before profile collision with ramp is %4.4f',EstFailTime);
    %     msgTimetoFailure
        %% 
        Xhybrid=Xprofile;
    %     if ( DXRight > 0) && (DXLeft > 0)
        %**********************************************************************
        %This is where the soliton is actually spliced on to shock-ramp profile
        %**********************************************************************
        if ProfileError == 0 || ProfileError == 2
            spliceSoliton = 1;
            if spliceSoliton ==1
                for ihybrid = 1:(idRsol - idLsol)
                    bKBW_IS(idxLsoliton + ihybrid-1) = -(Soliton.Locate)*profileSoliton(PrflSoliton, Physics,Xhybrid(ihybrid + idLsol-1),0);
                end
            end
            figure(BFieldFig);
            zmax = 1.5*max(abs(  bKBW_IS ) );
            zmin = -zmax;

            Proceed = 0;
            close all;

    %         bfit = bKBW_IS;           % bfit is default profile name to pass to boundary approval routine
            %Define figure for Hybrid IS - Soliton B-field Graph and Final Approval
            BFieldWindowPosition = [(NMonitors - 1 +0.1) 0.15 0.8 0.85]; 
            BFieldFigName = sprintf('%s Graph', PrflSoliton('Profile1'));  
            B1flag = 1;                                       % Changed to another value in while loop once user approves boundaries
            idL1 = 1;idR1 = N;                                % Assumption that profile is defined on entire range
            xL = -L;xR = +L;                                  % Assumption that profile is defined on entire range
            Xprofile = linspace(xL,xR,N);                     % Assumption that profile is defined on entire range 
            zmax = 1.5*max( abs(bfit) );                      % Determine y-scale for graphing
            zmin = -zmax;
            comment1 =' ';
            comment2 ='Proceed with this Hybrid Shock-Soliton Profile? ';
            FirstPass = 1;
            TxtAlign = 1;
            SolitonTxtAlign = +Soliton.Locate;
            ShockTxtAlign = -Soliton.Locate;


            while B1flag == 1

                BFieldFig = figure('Name',BFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BFieldWindowPosition,'visible','off');
                clf;

                axis([Xprofile(idL1) Xprofile(idR1) zmin zmax]);
                IS_Shock_Profile_Label(BFieldFig, Profile1, SolitonTxtAlign, zmax, xLsoliton, xRsoliton, Xprofile(idL1:idR1), bsol(idL1:idR1)   , comment1);
                IS_Shock_Profile_Label(BFieldFig, Profile2, ShockTxtAlign  , zmax, XShockL  , XShockR  , Xprofile(idL1:idR1), bKBW_IS(idL1:idR1), comment2);
                BFieldPlot(BFieldFig,bKBW_IS(idL1:idR1),Xprofile(idL1:idR1));
                [B1flag, FirstPass, Prfl_fit, xL_profile, xR_profile, L_profile, idL1, idR1] = Boundary_Approval(BFieldFig, Profile, Prfl_fit, 'Lsoliton', FirstPass, Xprofile, XShockL, XShockR);
                if B1flag ==1
                    close(BFieldFig);
                end

            end


                % Construct a questdlg with three options
                choice = questdlg('Proceed?', ...
                    'Proceed to numerical integration?', ...
                    'Yes','No', 'No');
                % Handle response
                switch choice
                    case 'Yes'
                        disp([choice 'Proceeding.'])
                        Proceed = 1;
                    case 'No'
                        disp([choice ' I''m out of here!'])
                        Proceed = 0;
                    case 'No'
                        disp([choice ' I''m out of here!'])
                end


        else
            fprintf(['\n\nNot possible with given value for L and given locations for the Shock and the Ramp.\n']);
            Proceed = 0;

    %         if CollisionDistance > 0
    %             fprintf(['\nIf Shock and Ramp locations are adjusted as indicated, run time might be as long as ',num2str(EstFailTime),' before failure.\n']);
    %             fprintf([msgIS,'\n']);
    %             fprintf([msgRmp,')L\n']);
    %         elseif CollisionDistance < 0
    %             msg = sprintf('Consider increasing L to about %4.4f for a valid run of est. time 200 -> 400\n',Ladvised);
    %             fprintf([msg,'\n']);
    %         end

        end

        zmax = 1.5*max(abs(bKBW_IS));
        zmin = -zmax;
        GraphicsValues = {Gtitle, exact, zmax,zmin, NMonitors};
        Graphics = containers.Map(GraphicsKeys, GraphicsValues);

        PrflSolitonKeys   = PrflSoliton.keys;
        PrflSolitonValues = PrflSoliton.values;

        PrflGeneratedISKeys   = PrflGeneratedIS.keys;
        PrflGeneratedISValues = PrflGeneratedIS.values;

        PrflKBWKeys   = PrflKBW.keys;
        PrflKBWValues = PrflKBW.values;

        PrflFinalKeys   = [PrflSolitonKeys  , PrflGeneratedISKeys  , PrflKBWKeys   ];
        PrflFinalValues = [PrflSolitonValues, PrflGeneratedISValues, PrflKBWValues ];
        Prfl = containers.Map(PrflFinalKeys , PrflFinalValues);
        PrflFinalKeystmp = Prfl.keys;
        PrflFinalValuestmp = Prfl.values;
        if Shock.Type == 1
            Prfl('Profile') =['1->3 IS-',Prfl('Profile1')];
        elseif Shock.Type == 2
            Prfl('Profile') =['2->3 IS lower branch-',Prfl('Profile1')];
        elseif Shock.Type == 3
            Prfl('Profile') =['2->3 IS upper branch-',Prfl('Profile1')];
        else
            Prfl('Profile') =['???',Prfl('Profile1')];
        end

        %                     DiagnosticsHybrdSplice = PStr.DCon;

        PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
            'PrflCon', Prfl, 'GCon',Graphics,'FCon',Folder, 'DrvrCon',Driver);
        Prfl = PStr.PrflCon;
        Prflkeys = Prfl.keys;
        Prflvals = Prfl.values;

        profile = @(Prfl, Physics,x,t) bKBW_IS;

    end
    %% ProfileType #12     Hybrid: Archived Trial vs. Soliton 
    if (ProfileType == 12)
        Prfltitle = 'Profile Parameters';
        Profile = 'Hybrid';
        Proceed = 0;
       %% Todo
    %   1) Load Parameters from original trial
    %   2) Be sure these parameters are saved in PStr ... and not overwritten
    %   ... apart from ...
    %   3) Update parameters in PStr to specify new:
    %       TotalTime, TimestepsAttempted, timestepsPerFrame, Courant ... others?
    %       while being sure that, L, N and physics remain the same.
    %   4) Be sure this PStr is passed to DNLS_v4
    %   5) Be sure original TrialID is saved to PStr

       %% Contruct Parameters for Archived Bfield Trial, Run and Line #'s
            TrialNumberOld = TrialNumber1; % Trial # of Archived B-field to use
            RunNumberOld   = RunNumber1;   % Run   # of Archived B-field to use
            LineNumberOld  = LineNumber1;  % Line  # of Archived B-field to use

    % Construct Parameters for New Bfield Trial and Run #'s
            TrialNumberNew = TrialNumber; % Trial # of B-field being created
            RunNumberNew   = RunNumber;   % Run   # of B-field being created

        %***********************************************************    

       %% Extract PStr and b1 from Archive of Old Trial
        %          ... As well as the specified B-field line, b1
        %***********************************************************     

    %     [PStr, b1, MSG, flag] = Read_Bfield_Line(ArchiveBase, TrialNumber1,RunNumber1, LineNumber1);
        TrialNumberTemp = 0;% set to 0 as a double, no Bfield lines are written to a new archive
                           % set to '17' (the trial # as a char string) ... the new trial number to populate with B-field lines prior to TrialNumberOld
        RunNumberTemp = 0;  % if TrialNumberNew = 0, the value for RunNumberNew is not used
        RBLdebug = 0;      % 0 --> no feedback is written to the command line
        [PStr, b1, LineNselected, MSG , flag ] = Continuation_Read_Bfield_Line(ActiveArchiveBase, TrialNumberOld, RunNumberOld, LineNumberOld, TrialNumberTemp, RunNumberTemp, RBLdebug);
        TrialHybrdSplice       = PStr.TCon;
        NumericsHybrdSplice    = PStr.NCon;
        PhysicsHybrdSplice     = PStr.PCon;
        PrflHybrdSplice        = PStr.PrflCon;
        GraphicsHybrdSplice    = PStr.GCon;
        FoldersHybrdSplice     = PStr.FCon;
        DriverHybrdSplice      = PStr.DrvrCon;
        DiagnosticsHybrdSplice = PStr.DCon;

       %% Update  Numerics Parameters: Folders Old and New, Courant, TotalTime, timestepsPerFrame, others?

       %***********************************************************************
       %  Update Numerics parameters as needed*

       L = NumericsHybrdSplice('L'); % Retrieved from Archived Parameters
       N = NumericsHybrdSplice('N'); % Retrieved from Archived Parameters
       dx = 2*L/(N-1);               % Recalculated based on Archived values
       dt = dx * Courant;            % May differ from Archived value if Courant is changed

       NumericsHybrdSplice('TotalTime') = TotalTime;
       NumericsHybrdSplice('Courant') = Courant;
       NumericsHybrdSplice('dx') = dx;
       NumericsHybrdSplice('dt') = dt;

       TimestepsAttempted = round(TotalTime/dt)+1; % Recalculated as TT and dt may differ from Archived values
       timestepsPerFrame  = round(TimestepsAttempted/NumberLinesSaved + .5); % Recalculated as TimestepsAttempted and NumberLinesSaved may differ

       NumericsHybrdSplice('TimestepsAttempted') = TimestepsAttempted; % Update given possible changes
       NumericsHybrdSplice('timestepsPerFrame')  = timestepsPerFrame;  % Update given possible changes

       %  Finished Updating Numerics parameters 
       %***********************************************************************

       %% Update Folder Parameters both Old and New
       %***********************************************************************
       %  Begin Updating Folder parameters: Old and New
       %***********************************************************************
       % Folder Paths For Archived Trial: TrialNumberOld & RunNumberOld

        TrialFolderRelOld = ['\Trial ',TrialNumberOld];%e.g. '\Trial 11'
        BfieldFolderRelOld = [TrialFolderRelOld,'\Bfields'];%e.g. '\Trial 11\Bfields'

        %Absolute Folder Paths
        TrialFolderAbsOld = [ActiveArchiveBase,TrialFolderRelOld];%e.g. 'F:\ ...\Trial 11\'
        BfieldFolderAbsOld = [ActiveArchiveBase,BfieldFolderRelOld];%e.g. 'F:\ ...\Trial 11\'

        %path to excel file to record trial parameters
        xlsfilenameOld = [NumericalMethodFile,'_',TrialNumberOld,'.xlsx'];% e.g. 'DNLS2D10c_11.xls'
        xlsfileOld = [TrialFolderAbsOld,'\',xlsfilenameOld];%full path to xls file

        %***********************************************************************
        % Folder Paths For New Trial: TrialNumberNew & RunNumberNew

        TrialFolderRelNew = ['\Trial ',TrialNumberNew];%e.g. '\Trial 11'
        BfieldFolderRelNew = [TrialFolderRelNew,'\Bfields'];%e.g. '\Trial 11\Bfields'

        %Absolute Folder Paths
        TrialFolderAbsNew = [ActiveArchiveBase,TrialFolderRelNew];%e.g. 'F:\ ...\Trial 11\'
        BfieldFolderAbsNew = [ActiveArchiveBase,BfieldFolderRelNew];%e.g. 'F:\ ...\Trial 11\'

        utilities(TrialFolderAbsNew, TrialFolderRelNew); %Creates Trial folder if needed
        utilities(BfieldFolderAbsNew, BfieldFolderRelNew); %Creates Bfield folder if needed

        %path to excel file to record trial parameters
        xlsfilenameNew = [NumericalMethodFile,'_',TrialNumberNew,'.xlsx'];% e.g. 'DNLS2D10c_11.xls'
        xlsfileNew = [TrialFolderAbsNew,'\',xlsfilenameNew];%full path to xls file

        MatID   = strcat(TrialNumberNew);
        B_field_Trial_Run_matfile_rel = ['\Bfield_',MatID,'.mat'];%
        B_field_Tr_mat     = [BfieldFolderAbsNew, B_field_Trial_Run_matfile_rel];

        TrialIDNew   = strcat(TrialNumberNew,'_',int2str(RunNumberNew)); % "#-#' format
        TrialK_B_field_relNew = ['\Bfield_',TrialIDNew,'.txt'];
        BfieldfileNew      = [BfieldFolderAbsNew, TrialK_B_field_relNew];    

        % Store KeysOld to store OLD Folder values
        FoldersHybrdSplice('BfieldFolderOld')  = FoldersHybrdSplice('BfieldFolder');
        FoldersHybrdSplice('BfieldMatfileOld') = FoldersHybrdSplice('BfieldMatfile');
        FoldersHybrdSplice('BfieldfileOld')    = FoldersHybrdSplice('Bfieldfile');
        FoldersHybrdSplice('TrialFolderOld')   = FoldersHybrdSplice('TrialFolder');

        % Store NEW Folder values as default Keys
        FoldersHybrdSplice('BfieldFolder')  = BfieldFolderAbsNew;
        FoldersHybrdSplice('BfieldMatfile') = B_field_Tr_mat;
        FoldersHybrdSplice('Bfieldfile')    = BfieldfileNew;
        FoldersHybrdSplice('TrialFolder')    = TrialFolderAbsNew;    

       %  Finished Updating New and Old Folder Paths 
       %***********************************************************************

       %% Update Trial Parameters: Shift NEW to 'Trial' and OLD to 'TrialOld'

       % Store Old Parameters
       TrialHybrdSplice('RunNumberOld')   = RunNumberOld;
       TrialHybrdSplice('TimeStampOld')   = TrialHybrdSplice('TimeStamp');
       TrialHybrdSplice('TrialIDOld')     = TrialHybrdSplice('TrialID');   
       TrialHybrdSplice('TrialNameOld')   = TrialHybrdSplice('TrialName');
       TrialHybrdSplice('TrialNumberOld') = TrialNumberOld;

       % Create Parameters for current trial
       TrialIDNew   = strcat(TrialNumberNew,'_',int2str(RunNumberNew)); % "#-#' format
       TrialNameNew = [NumericalMethodFile,'_',TrialIDNew]; % 'DNLS2D10c_11_4'
       TimeStamp = datestr(clock,0);
       TrialK_B_field_relNew = ['\Bfield_',TrialIDNew,'.txt'];

       % Save Parameters for current trial
       TrialHybrdSplice('RunNumber')   = int2str(RunNumberNew);
       TrialHybrdSplice('TimeStamp')   = TimeStamp;
       TrialHybrdSplice('TrialID')     = TrialIDNew;   
       TrialHybrdSplice('TrialName')   = TrialNameNew;
       TrialHybrdSplice('TrialNumber') = TrialNumberNew;
       TrialHybrdSplice('TrialObjective') = 'Collide a non coplanar intermediate shock with a soliton.';

       %% Interpolate Bfield, B1, of Shock to match original grid
        X4uniform    = linspace(-L, +L, length(b1));
        Xinterpolate = linspace(-L, +L, N);   
        B1 = interp1(X4uniform, b1, Xinterpolate); % B-field from old trial to be passed to DNLS_v4 as start of extended trial
        X1 = Xinterpolate;

       %% Identify Shock Boundaries with First and Second Click and with Third Click identify Quiescent Field on Either Left or Right Side of Shock Where Soliton is to be Spliced
        %*************************************************************************************************************    
        %      Graph B1: Archived B-field to splice soliton into
        %      Select: --> the left & right asymptotes of the Archived B-field
        %              --> Last, select the background field on side into which
        %              the soliton will be spliced.
        %
        %                                     |   Shock Section   | 
        %                                     |         xx        |
        %                                     |      xxx  x       |
        %                                     |   xxx      x      |
        % The FIRST Click at LHS of Shock --> |xxx         x   xx |<-- The SECOND Click at RHS of Shock
        %                                     |             x  x  |
        %                                     |              xx   |
        %
        %   The THIRD and last Click: 
        %   --> Click LEFT  of Shock in a spot where the B-field is constant to place the soliton on the LEFT  side
        %   --> Click RIGHT of Shock in a spot where the B-field is constant to place the soliton on the RIGHT side          
        %*************************************************************************************************************    
        close all;
        B1WindowPosition = [(NMonitors - 1 +0.1) 0.1 0.8 0.8];
        B1FieldFigName = sprintf('Shock-Soliton Configuration');
        B1FieldFig = figure('Name',B1FieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',B1WindowPosition,'visible','off');

        plot(X1,real(B1), 'g');
        hold('on');
        plot(X1,imag(B1), 'r');
        hold('on');
        plot(X1,abs(B1), 'b');
        hold('on');

        profilename = PStr.PrflCon('Profile'); % Returns profile name as a cell array, but PrflCon('Profile') returns a char array
        text(-L*0.9,0.92*zmax,sprintf('Profile = %s',profilename));
        text(-L*0.4,0.92*zmax,sprintf('Trial ID = %s',PStr.TCon('TrialID')));
        text(-L*0.9,0.85*zmax,sprintf('L = %4.4f',PStr.NCon('L')  ));
        text(-L*0.4,0.85*zmax,sprintf('N = %4.4f',PStr.NCon('N')  ));
        text(-L*0.9,0.78*zmax,sprintf('line # = %4.0f',LineNumber1));

        axis([-L L zmin zmax]);
        drawnow();

        B1FieldFig.Visible = 'on';
        B1flag = 1;
        while B1flag == 1
            [x,y] = ginput(3);
            xL = x(1);
            xR = x(2);
            xAS = x(3);

            if xR < xL      % Switching xL and xR if reversed
                xtemp = xR;
                xR = xL;
                xL = xR;
            end

            SolitonLocation = +0;
            if xAS > xR  % Soliton on RHS of Shock
                SolitonLocation = +1;
            else
                SolitonLocation = -1;
            end

            B1FieldFig;

            plot([xAS xAS],[zmin zmax]*0.8,'Color',[0 0 1]);
            hold('on');
            text(xAS + 2*L*0.01,0.8*zmin,'Soliton inserted on this side?');

            tmpL = abs(X1 - xL); tmpR = abs(X1 - xR); tmpAS = abs(X1 - xAS);
            [val, idLnew] = min(tmpL); [val, idRnew] = min(tmpR); [val, idAS] = min(tmpAS);
            if xAS > xR    % Soliton and xAS to the RHS
                bR1 = B1(idAS);   % B-field on the RHS
                bL1 = B1(idLnew); % B-field on the LHS
            else            % Soliton and xAS to the LHS
                bR1 = B1(idRnew); % B-field on the RHS
                bL1 = B1(idAS);   % B-field on the LHS
            end

            %Evaluate B-Field Rotation
            e1R = [real(bR1) imag(bR1)]/abs(bR1);
            e1L = [real(bL1) imag(bL1)]/abs(bL1);
            angle = acos(e1L*e1R')*180/3.14159265358979;
            comment1 = sprintf('Shock and Soliton Locations: Field Rotation %4.2f degrees.',angle);
            IS_Shock_Profile_Label(B1FieldFig, 'Shock', 0, zmax, xL, xR, X1, B1, comment1) ;

            % Construct a questdlg with three options
            choice = questdlg('Selection correct?', ...
                'Asymptotes for profile and transition to background state:', ...
                'Yes','Try Again','Get Me Outa Here','Get Me Outa Here');
            % Handle response
            switch choice
                case 'Try Again'
                disp([choice ' Try, try again.']);
                B1flag = 1;
                case 'Yes'
                disp([choice ' Woot!'])
                B1flag = 2;
                case 'Get Me Outa Here'
                disp('Prepare for Transport.')
                B1flag = -1;
            end
        end


        %***********************************************************    
        % Create B-field lines for B2 (taken to be a 1p-Soliton)
        % Define Prfl2 for lambda, sign and bl defining the 1p-Soliton
        %***********************************************************
        b1p = B1(idAS);
        b1 = b1p;

       %% Define Prfl2 for Soliton and join with PrflHybrdSplice
        lambda = Soliton.lambda;   % Default lambda set here s.t. 0 < lambda < 1 
        sign1 = Soliton.sign;       % sign = +1 ... 'bright soliton'; sign = -1 ... for a 'dark soliton'
        Profile2 = '1psoliton';

        %*********************************************************
        v = b1^2 + 2*lambda^2;  % the soliton speed ... to be updated with changes above
        PrflKeys2   = {'title2'   , 'Profile2' , 'speed2', 'lambda',  'sign', 'b1p'};
        PrflValues2 = { Prfltitle,  Profile2 ,    v   ,  lambda ,   sign1 ,  b1p };
        Prfl2       = containers.Map(PrflKeys2,PrflValues2);
        %*********************************************************  

        %*********************************************************
        %    Joing Shock and Soliton Prfl's
        PrflHybrdSplice('Profile') = 'KBW-Soliton Spliced';
        PrflHybrdSpliceKeys2   = [PrflHybrdSplice.keys, Prfl2.keys];
        PrflHybrdSpliceValues2 = [PrflHybrdSplice.values, Prfl2.values];
        PrflHybrdSplice = containers.Map(PrflHybrdSpliceKeys2, PrflHybrdSpliceValues2); 

        %*********************************************************

       %% Identify Soliton Boundaries
        B2 = onePsoliton(Prfl2, PhysicsHybrdSplice, X1, 0);    
         B2asymptote = abs(B2(1) - b1p) + abs(B2(end) - b1p);
         if B2asymptote > 0.02
             B2FieldFigName = sprintf('The soliton may be too wide for the width 2L. Try a larger lambda?');
         else
             B2FieldFigName = sprintf('Soliton Configuration');
         end
        close(B1FieldFig);
        B2WindowPosition = [(NMonitors - 1 +0.1) 0.1 0.8 0.8];

        B2FieldFig = figure('Name',B2FieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',B2WindowPosition,'visible','off');

        B2FieldFig;
        plot(X1,real(B2), 'g');
        hold('on');
        plot(X1,imag(B2), 'r');
        hold('on');
        plot(X1,abs(B2), 'b');
        hold('on');

        profilename = Prfl2('Profile2');
        text(-L*0.9,0.9*zmax,sprintf('Profile = %s',profilename));
        text(-L*0.4,0.89*zmax,sprintf('bl = %4.4f',Prfl2('b1p')));
        text(-L*0.9,0.79*zmax,sprintf('lambda = %4.4f',Prfl2('lambda')  ));
        text(-L*0.4,0.79*zmax,sprintf('sign = %1.0f',Prfl2('sign')  ));

        axis([-L L zmin zmax]);
        drawnow();

        B2flag = 1;
        B2FieldFig.Visible = 'on';
        while B2flag == 1
            [x,y] = ginput(2);
            x = sort(x);
            x2L = x(1);x2R = x(2);
                  tmp2L  = abs(X1 - x2L);      tmp2R  = abs(X1 - x2R)   ;
            [val, idx2L] = min(tmp2L)   ;[val, idx2R] = min(tmp2R);

            bsolAsymptote = (B2(idx2L) + B2(idx2R))/2; % Averaging assuming selected locations may vary slightly from background

            B2FieldFig;
            %Evaluate Soliton Orientation
            e2R = [real(bsolAsymptote) imag(bsolAsymptote)]/abs(bsolAsymptote);
            e2L = [1 0];
            angle2 = acos(e2L*e2R')*180/3.14159265358979;
            B2asymptote = abs(B2(1) - b1p) + abs(B2(end) - b1p);
            if B2asymptote > 0.02
                comment2 = sprintf('The soliton may be too wide for the width 2L. Try a larger lambda?');
            else
                comment2 = sprintf('Locating Soliton Boundaries: Soliton oriented %4.2f degrees from vertical',angle2);
            end

            IS_Shock_Profile_Label(B2FieldFig, 'Shock', 0, zmax, x2L, x2R, X1, B2, comment2) ;

            % Construct a questdlg with three options
            choice = questdlg('Selection correct?', ...
                'Select Section of Profile B2 to use:', ...
                'Yes','Try Again','Get Me Outa Here','Get Me Outa Here');
            % Handle response
            switch choice
                case 'Try Again'
                disp([choice ' Try, try again.'])
                B2flag = 1;
                case 'Yes'
                disp([choice ' Woot!'])
                B2flag = 2;
                case 'Get Me Outa Here'
                disp('Prepare for Transport.')
                B2flag = -1;
                Proceed = 0;
            end
        end    

       %% Construct Hybrid Profile
        data = B1;  % 'data' will be the final hybrid profile and is initialized to be the original shock profile
        % Determine if Soliton was placed on the LHS or the RHS of shock
        if     SolitonLocation == +1
            INLeft  = idRnew;
            INRight = idRnew + (idx2R - idx2L);
        elseif SolitonLocation == -1
            INLeft  = idLnew - (idx2R - idx2L);
            INRight = idLnew;
        else
            ['Something went sideways in locating the soliton.']
            return;
        end

        j2 = idx2L;
        for j = INLeft:INRight
            data(j) = B2(j2);
            j2 = j2 + 1;
        end  % This splices the soliton into the selected section of the shock profile

        close(B2FieldFig);
        BHWindowPosition = [(NMonitors - 1 +0.1) 0.1 0.8 0.8];
        BHFieldFigName = sprintf('Final Profile');
        BHFieldFig = figure('Name',BHFieldFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',BHWindowPosition,'visible','off');

        BHFieldFig;
        plot(X1,real(data), 'g');
        hold('on');
        plot(X1,imag(data), 'r');
        hold('on');
        plot(X1,abs(data), 'b');
        hold('on');

        profilename = 'Hybrid';
        text(-L*0.9,0.9*zmax,sprintf('Profile = %s',profilename));
        text(-L*0.4,0.89*zmax,sprintf('bl = %4.4f',Prfl2('b1p')));
        text(-L*0.9,0.79*zmax,sprintf('lambda = %4.4f',Prfl2('lambda')  ));
        text(-L*0.4,0.79*zmax,sprintf('sign = %1.0f',Prfl2('sign')  ));

        axis([-L L zmin zmax]);
        drawnow();

        Hybrid.Select = 0;
        BHFieldFig.Visible = 'on';
        % Construct a questdlg with three options
        choice = questdlg('Proceed with this Profile?', ...
            'Profile Approval:', ...
            'Yes','Try Again','Get Me Outa Here','Get Me Outa Here');
        % Handle response
        switch choice
            case 'Yes'
                disp([choice ' To infinity and beyond!'])
                Hybrid.Select = 1;
                Proceed = 1;
            case 'Try Again'
                disp([choice ' Not handled yet.'])
                Hybrid.Select = -1;
            case 'Get Me Outa Here'
                disp('Not handled yet.')
                Hybrid.Select = -2;
                Proceed = 0;
        end
        close all;
        profile = @(Prfl, PhysicsHybrdSplice,x,t) data;

    end

    %% *******************************************
    %% Execute: Numerical Method & Archiving Features
    %********************************************************************
    %% Verify that Trial does not already exist

    if (exist(B_field_Tr_Rn_txt,'file') == 2 )  && Proceed
        BFtxtparts = split(B_field_Tr_Rn_txt,{'\'});
        fileIdx = length(BFtxtparts);
        Question = sprintf('%s already exists.',BFtxtparts{fileIdx});

        choice = questdlg(Question, ...
            'Overwrite?                                                         ', ...
            'Yes','No','No');
        % Handle response
        switch choice

            case 'Yes'
                Proceed = 1;

            case 'No'
                Proceed = 0;
        end
    end
    %% B-Field Archive Header
    if (BFieldarchive == 1) && Proceed == 1
        %Write B-field Archive File Header

        elapsed = num2str(cputime - t); t = cputime;
        ['Starting to write B-field Header: ', TrialName, ' --- Elapsed time = ',elapsed, ' s']
        Prfl = PStr.PrflCon;
        Prflkeys = Prfl.keys;
        Prflvals = Prfl.values;
        DNLS_Bfield_Header_v2(PStr,'w', RunNumber, t);
        Prfl = PStr.PrflCon;
        Prflkeys = Prfl.keys;
        Prflvals = Prfl.values;
        elapsed = num2str(cputime - t); t = cputime;
        ['Completed writing B-field Header: ', TrialName, ' --- Elapsed time = ',elapsed, ' s']
    end
    %% Numerical Method
    if (runNumericalMethod == 1) &&( Proceed ==1)
        %Call to Numerical Method
        elapsed = num2str(cputime - t); t = cputime;
        CPUtime = num2str(0.00) ;timesteps = 1; growth = 0; decay = 0.9;
        ['Starting: ', TrialName, ' --- Elapsed time = ',elapsed, ' s']
        Prfl = PStr.PrflCon;
        Prflkeys = Prfl.keys;
        Prflvals = Prfl.values;

        [PStr, growth, decay, timesteps] = DNLSHandle(profile, PStr);
        Prfl = PStr.PrflCon;
        Prflkeys = Prfl.keys;
        Prflvals = Prfl.values;

        CPUtime = num2str(cputime - t); t = cputime;
        ['Completed:', TrialName, ' --- Elapsed time = ',CPUtime, ' s']
    end
    %% Archiving
    %Prepare Post Trial Diagnostic Parameters
    if (Proceed == 1)
        if (runNumericalMethod ~= 1)    
            CPUtime = num2str(0.00) ;timesteps = 1; growth = 0; decay = 0.9;
        end
        Dtitle = 'Post Run Diagnostics';
        DiagnosticValues = {Dtitle,CPUtime, timesteps, growth, decay};
        DiagnosticKeys = {'title', 'CpuTimeUsed', 'timesteps','growth', 'decay'};
        Diagnostics = containers.Map(DiagnosticKeys, DiagnosticValues);
        Prfl = PStr.PrflCon;
        PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
            'PrflCon', Prfl,'GCon',Graphics,'FCon',Folder,'DrvrCon',Driver,'DCon',Diagnostics);
    end
    %**********************************************************

    if ( Excelarchive == 1 ) && ( Proceed == 1 )
        %Write Trial Parameters to Excel Archive
        Prfl = PStr.PrflCon;
        Prflkeys = Prfl.keys;
        Prflvals = Prfl.values;
        DNLS_Excel_Archive_v2(xlsfile, PStr, RunNumber, t)
        %             DNLS_Excel_Archive(xlsfile, PStr, k, t)
    end
    if (BFieldarchive == 1) &&( Proceed ==1)
        % Create BFfileMat
        if CreateMatfile == 1
            BFMat = matfile(B_field_Tr_mat,'Writable',true);
            BFMat.PStr(RunNumber,1) = PStr;
            clear BFMat;
        end
        %Write B-field Archive File Footer
        elapsed = num2str(cputime - t); t = cputime;
        ['Starting to write B-field Footer: ', TrialName, ' --- Elapsed time = ',elapsed, ' s']
        DNLS_Bfield_Header_v2(PStr, 'a', RunNumber, t)
        elapsed = num2str(cputime - t); t = cputime;
        ['Completed writing B-field Footer: ', TrialName, ' --- Elapsed time = ',elapsed, ' s']
    end
    if (SendEmail == 1) &&( Proceed ==1) && (  RunNumber == RunNumber_end )
        %Send email with xlsfile attached
        loginfile = strcat(pwd,'\Bobs_email_data.xls');
        elapsed = num2str(cputime - t); t = cputime;
        ['Sending email. ','Elapsed time = ',elapsed, ' s']
        sendmailguru(loginfile,xlsfile, TrialID, TimeStamp)
        elapsed = num2str(cputime - t); t = cputime;
        ['Email Sent. ','Elapsed time = ',elapsed, ' s']
    end
end





if PlaySound ==1
    load handel;
    sound(y,Fs);
end
%% Function Definitions
    function [XrampNew, XISNew, TBCFail, Tcollision] = HybridOptimize(Lo, LPrflo, Lrampo, Lsolitono, LShocko, Sol, vsolo, vis, vrw, VRWFooto)
        DXfree = 2*Lo - LPrflo;
        impactSpeedo = vis - vsolo;
        Tcollision = 1.1*(LShocko + Lsolitono)/impactSpeedo;
        
        % numerically it is observed that the expression used above for vrw
        % is significantly too slow. 
        % With bL = 1 and bk = 0.25 it was observed, roughly, 
        % that vrw ~ 3*vis
        
        
        % The free distance figured above, DXfree, needs to be split between ...
        % --> dxr = buffer distance between the right of the shock-soliton and the ramp 
        % --> dxl = buffer between the left of the shock-soliton and the right edge of the ramp.
        %
        % Note dxr != dxl as dxr = (vis - VRWFooto)*time and dxl = vrw*time and the simulation time is the same for both
        % Therefore DXfree = dxr + dxl = (vis - VRWFooto + vrw)*time
        % --> TBCFail = DXfree/(vis + vrw) ... is the time at which the
        % physics boundary conditions fail
        % dxr = vis*DXfree/(vis - VRWFooto + vrw)
        % dxl = vrw*DXfree/(vis - VRWFooto + vrw)
        %
        %  left edge of simulation --> | ... dxl/2 ...| Shock-Soliton| ...dxr ...| Ramp| ... dxl/2 ... | <-- right edge of simulation

        TBCFail = DXfree/(vis - VRWFooto + vrw);
        dxr = (vis - VRWFooto)*DXfree/(vis - VRWFooto + vrw);
        dxl = vrw*DXfree/(vis - VRWFooto + vrw);
        XrampRnew = Lo - dxl/2;
        XrampNew = XrampRnew - Lrampo/2;
        XrampLnew = XrampNew - Lrampo/2;

        if Sol.Locate == -1       % soliton to the left and shock to the right
            XISRight = XrampLnew - dxr;
            XISNew = XISRight - LShocko/2; % Intermediate shock has been observed to be roughly centered at specified XIS

        elseif Sol.Locate == +1   % shock to the left and soliton to the right
            XSSRight = XrampLnew - dxr;
            XISNew = XSSRight - Lsolitono - LShocko/2; % Intermediate shock has been observed to be roughly centered at specified XIS

        end
    end
    
    function [out] = Si(x) %Pade approximation for Sine Integral (see Appendix eq. (B.5) of Rowe et. al.: https://arxiv.org/pdf/1407.7676.pdf)
        
        num = 1 - 4.54393409816329991*10^(-2)*x^2 + ...
        1.15457225751016682*10^(-3)*x^4 - 1.41018536821330254*10^(-5)*x^6 + ...
        9.43280809438713025*10^(-8)*x^8 - 3.53201978997168357*10^(-10)*x^10 + ...
        7.08240282274875911*10^(-13)*x^12 - 6.05338212010422477*10^(-16)*x^14;
        
        den = 1 + 1.01162145739225565*10^(-2)*x^2 + ...
        4.99175116169755106*10^(-5)*x^4 + 1.55654986308745614*10^(-7)*x^6 + ...
        3.28067571055789734*10^(-10)*x^8 + 4.5049097575386581*10^(-13)*x^10 + 3.21107051193712168*10^(-16)*x^12;
        
        out = x*num/den;
    end
    
    function BFieldPlot(BFig,bb,x)
        figure(BFig);
%         clf;
        Bmax = 1.5*max(abs(bb));
        Bmin = -Bmax;
        plot(x,real(bb),'Color',[ 0, 1, 0]);hold 'on'; 
        plot(x,imag(bb),'Color',[ 1, 0, 0]); hold 'on';
        plot(x,abs(bb),'Color',[ 0, 0, 1]);hold 'off';
        axis([x(1), x(end), Bmin, Bmax]);
        drawnow();
    end
    
    function[idMax  , idSL     , idSR     , xShockL, xShockR, Iz    ] = IS_ShockBoundaries(Profile, Sol    , Shck , bb  , x       , bL         , basymptoteL, basymptoteR, ntest       , tol        )
    %%     [idRelMax, idxShockL, idxShockR, XShockL, XShockR, Izcalc] = IS_ShockBoundaries(Profile, Soliton, Shock, bfit, Xprofile, b_soliton  , 0          , 0          , NTestSoliton, Tol_Soliton);
    
    %% ARGUMENTS: bb, x, bL, basymptote, ntest, tol
        % HybIS.Profile  '0'--> No soliton to be spliced; '1' --> Soliton to be spliced
        % HybIS.Locate  '-1'--> place soliton to the left of shock; '+1' --> place soliton to the right of the shock
        % HybIS.Soliton  '1'--> 1p soliton; '2' --> 2p soliton
        % HybIS.Shock  '0'--> No shock; '1' --> 13 IS, '2' --> 23 IS with lower Iz, '3' --> 23 IS with higher Iz
        
%         Soliton.Type  = 1   ;  % 1 --> 1p soliton, 2 --> 2p soliton, 3 --> KNS soliton
%         Soliton.lambda = 0.02;  % For a 1p soliton & actual eigenvalue is lambda_actual = lambda*bo
%         Soliton.sign  = 1  ;  % +1 --> bright 1p soliton, -1 --> dark 1p soliton
%         Soliton.lr    = 0.5 ;  % Real part of eigenvalue for a 2p soliton --> actual eigenvalue is lambda_actual = (lr + I*li) *bo
%         Soliton.li    = 0.4 ;  % Imaginary part of eigenvalue for a 1p soliton
%         Soliton.Locate = -1; % -1 --> place soliton on LHS of shock, +1 --> place soliton on RHS of shock

%         Shock.Type  =  1    ;  % 1 --> 13 IS, 2 --> 23 IS on lower branch, 3 --> 23 IS on upper branch
%         Shock.CR    = -0.42 ;  % Compression ratio of 13 IS with -1 < Shock.CR < 0 and CR23 = -1 - IS.CR13
%         % For Fast Shock 0 < CR < 1
%         Shock.Iz    = - 1.52; % Target Iz for IS
% 
%         % theta_L and theta_kbw are only used in the KBW profile
%         Shock.theta_L = 0*pi; % Normally set at 0 to match a Bfield in the vertical direction at the far LHS (by > 0 and bz = 0)
%         Shock.theta_kbw = 175/180*pi*(- sign(Shock.Iz) ); % Direction of Bfield downstream of shock (only used in KBW)
% 
%         % Specify location for shock and ramp
%         Shock.x     = -0.4395*L;  % location of midpoint of Shock in range: -L --> +L
%         Shock.xRamp = +0.3469*L;  % location of ramp ... in range: -L --> +L
        
        % bb --> section of B-field with IntermediatE Shock on LHS preceeded
        % by section with constant B-field = ( bL, 0)
        %
        % x --> array of position values of same length as bb
        %
        % bL ---> as described in bb above
        %
        % basymptote --> it is assumed that the z-component of the B-field
        % goes to this constant asymptotic value to the right of the shock.
        % It is should be evaluated prior to call as basymptote = b_L*sin(theta_kbw)
        %
        %ntest --> number of indices over which the running average of Bz
        %is evaluated   
        
    %% OUTPUTS: idMax, idSL, idSR, xShockL, xShockR, Iz
        % idMax --> the location of the maximum value of the z-component of
        % the B-field. For a dispersive Intermediate Shock the z-component
        % eventually settles into an oscillatory wave similar in appearance
        % to a damped sine-wave. prior to settling to this steady-state
        % this is not a well defined measure of shock's location.
        %
        % idSL --> index of the LHS of Shock
        % idSR --> index of the RHS of Shock
        % xShockL --> postion of the LHS of Shock
        % xShockR --> position of the RHS of Shock
        
        %evaluate the index of a 'fixed' location within the shock
        %Note: this location is only fixed once the shock reaches
        %steady-state
        
    %% Find first peak of Im(bb): idMax starting point for search of RHS
    if contains(Profile, 'Soliton')
        if Sol.sign == -1
            [val, idMax] = min( abs( bb ) ); % for dark 1p soliton, locates central minimum
        elseif Sol.sign == 1
            [val, idMax] = max( abs( bb ) ); % for bright 1p soliton, locates central maximum
        end
    else
        [val, idMax] = max(abs(imag(bb)));   % Assumed that alternative is an IS
    end
        
        
    %% Find LHS of Shock: 
    %     1st) assume Im(b) has monotonic decrease to asymptote to the left
    %     2nd) see if Im(b) has oscillatory approach to asymptote
        % Check for monotonic decrease to asymptote on LHS
        
        BLasymptoteSuccess = 0;
        minDeltaBz = 10;       
        idSL = 1;
        Lexception = 0;
        for idb = 1:(idMax - (ntest - 1))
            idx = idMax-(idb-1);
%             DeltaBz = abs(imag(bb(idx)) - basymptoteL);
%             
%             if DeltaBz < minDeltaBz
%                 minDeltaBz = DeltaBz;
%                 Lexception = 0;
%             else
%                 Lexception = Lexception +1';
%             end
%             
%             if Lexception > 50
%                 msg = 'Not monotonic'
%                 minDeltaBz = 10;
%             end
            
           
            idxLeft = idx-(ntest - 1);
            idxRight = idx + (ntest - 1);
            
            DeltaBzAvgL = sum(imag(bb(idxLeft:idx)))/ntest - basymptoteL;
            DeltaBzAvgR = sum(imag(bb(idx:idxRight)))/ntest - basymptoteL;
            
            SGNBzAvgL = sign(DeltaBzAvgL);
            SGNBzAvgR = sign(DeltaBzAvgR);
            
            bbzslopeR = (imag(bb(idxRight)) - imag(bb(idx)))/(x(idxRight) - x(idx));
            bbzslopeL = (imag(bb(idx)) - imag(bb(idxLeft)))/(x(idx) - x(idxLeft));
            
            SGNCrossing = SGNBzAvgL*SGNBzAvgR; % both should have the same sign if not evaluated at a crossing
            SGNslopes = sign(bbzslopeR*bbzslopeL);% Should be positive if approaching an asymptote from either above or below.
                                                  % Will be negative if coincidentally at min or max
            if (abs(DeltaBzAvgL) < tol*abs(bL)) &&  (abs(DeltaBzAvgR) < tol*abs(bL))  && (SGNCrossing > 0) && (SGNslopes > 0)
                idSL = idx;
                BLasymptoteSuccess = 1;
                break;
            end
        end
        
        BLoscillatory = 0;
        if ~BLasymptoteSuccess
            idSL = 1;
            for idb = 1:(idMax - ntest + 1)
                idxLeft = idMax-(idb-1)-(ntest - 1);
                idxRight = idMax-(idb-1);
                BzAvg = sum(imag(bb(idxLeft:idxRight)))/ntest;
                signLhsidb = sign(imag(bb(idxLeft)) - basymptoteL);
                signLhsidbPlusone = sign(abs(bb(idxLeft+1)) - basymptoteL);
                if (( signLhsidb*signLhsidbPlusone ) < 0) && (abs(BzAvg - basymptoteL ) < tol*abs(bL))
                    idSL = idb;
                    BLoscillatory = 1;
                    break;
                end
            end
        elseif ~BLoscillatory && ~BLasymptoteSuccess
            idSL = floor(idMax/2);
            msg = 'No well-defined LHS for Shock: ... going with midpoint'
        end
        
        xShockL = x(idSL);
    
    %% Find RHS of Shock: assumes Im(b) has general form of damped sign wave towards asymptote
        %Assuming Coplanar shock and imag(Bh) --> 0 to right of shock
        Nb = length(bb);
        idSR = floor((Nb + idMax)/2);
        Rexception = 0;
        
        % Check for monotonic decrease to asymptote on LHS
        BRasymptoteSuccess = 0;
        minDeltaBz = 10;
        for idb = (idMax + (ntest - 1)):(Nb - (ntest - 1))
%             DeltaBz = abs(imag(bb(idb)) - basymptoteR);
%             if DeltaBz < minDeltaBz
%                 minDeltaBz = DeltaBz;
%                 Rexception = 0;
%             else
%                 Rexception = Rexception +1;
%             end
%             
%             if Rexception > 50
%                 msg = 'Not monotonic'
% %                 break;
%             end
            
            idxLeft = idb - (ntest - 1);
            idxRight = idb + (ntest - 1);
             
            DeltaBzAvgL = sum(imag(bb(idxLeft:idb)))/ntest - basymptoteR;
            DeltaBzAvgR = sum(imag(bb(idb:idxRight)))/ntest - basymptoteR;
            
            SGNBzAvgL = sign(DeltaBzAvgL);
            SGNBzAvgR = sign(DeltaBzAvgR);
            
            bbzslopeR = (imag(bb(idxRight)) - imag(bb(idb)))/(x(idxRight) - x(idb));
            bbzslopeL = (imag(bb(idb)) - imag(bb(idxLeft)))/(x(idb) - x(idxLeft));
            
            SGNCrossing = SGNBzAvgL*SGNBzAvgR; % both should have the same sign if not evaluated at a crossing
            SGNslopes = sign(bbzslopeR*bbzslopeL);% Should be positive if approaching an asymptote from either above or below.
                                                  % Will be negative if coincidentally at min or max
           
            if (abs(DeltaBzAvgL) < tol*abs(bL)) &&  (abs(DeltaBzAvgR) < tol*abs(bL))  && (SGNCrossing > 0) && (SGNslopes > 0)
                idSR = idb;
                BRasymptoteSuccess = 1;
                break;
            end
        end

        BRoscillatory = 0;        
        if ~BRasymptoteSuccess
            for idb = idMax:(Nb - (ntest - 1) )
                idxLeft = idb;
                idxRight = idb + (ntest -1);
                Bzavg = sum( imag( bb( idxLeft:idxRight ) ) )/ntest;
                if ( ( imag(bb(idb-1)) - basymptoteR )*( imag(bb(idb)) - basymptoteR )) < 0 && (abs(Bzavg - basymptoteR ) < tol*abs(bL))
                    idSR = idb;
                    BRoscillatory = 1;
                    break;
                end
            end
        end
        if ~BRoscillatory && ~BRasymptoteSuccess
                idSR = floor(Nb);
                msg = sprintf('No well defined RHS for %s: ... going with far RHS./n Probably should zoom out.',Profile);
                msg
        end
        xShockR = x(idSR);
        
        Bis = bb(idSL:idSR);
        Xis = x(idSL:idSR);
        dxx = (Xis(1:end-1) - Xis(2:end)); % dxx defined as '-' to be consistent with convention in KBW paper in integrating across the shock from right to left
        Iz = dot(imag(Bis(2:end)),dxx);
        
    end
    
    function IS_Shock_Profile_Label(FigHandle, Profile, HFlag, Zmax, xSL, xSR, x, b, comment)
        
        % Invoke Figure
        figure(FigHandle);
        Nmax = length(b);
%         Zmax = 1.5*max(abs(  b ) );
        Zmin = -Zmax;
        idOS = floor(0.05*Nmax);
        tmpL = abs(x - xSL); tmpR = abs(x - xSR);
        [val, idSL] = min(tmpL); [val, idSR] = min(tmpR); 
        LShock = xSR - xSL;
        
        idCommTxtR = floor((idSL + idSR)/2);              % Place 'comment' at upper LHS
        XCommTxt = x(idCommTxtR); YCommTxt = 1.05*Zmax;
        text(XCommTxt,YCommTxt,sprintf('%s ',comment),'Color',[0 0 0],'horizontalAlignment', 'center', 'FontSize', 18);
        hold 'on';
        
        
        axis([x(1), x(end), Zmin, Zmax]);
        pos = get(gca, 'Position');
        NormShockXL=(xSL+abs(min(xlim)))/diff(xlim)*pos(3) + pos(1);% +AX(1)/AXrange;
        NormShockXR=(xSR+abs(min(xlim)))/diff(xlim)*pos(3) + pos(1);
        NormShockY=(0.5*Zmin - Zmin)/(Zmax - Zmin)*pos(4) + pos(2);% +AX(3)/AYrange;
        NormShockYTxt=(0.75*Zmin - Zmin)/(Zmax - Zmin)*pos(4) + pos(2);% +AX(3)/AYrange;
%         if contains(Profile,'Shock')
%             NormShockYTxt=(0.75*Zmin - Zmin)/(Zmax - Zmin)*pos(4) + pos(2);% +AX(3)/AYrange;
%         else
%             NormShockYTxt=(0.95*Zmin - Zmin)/(Zmax - Zmin)*pos(4) + pos(2);% +AX(3)/AYrange;
%         end
        NormShockYTxtHgt=0.10;% +AX(3)/AYrange;
        NormShockYTxtWdth=(xSR - xSL)/diff(xlim)*pos(3) ;% +AX(3)/AYrange;
        
        alignment = 'center';
        NormShockTxtXL = NormShockXL;
        if HFlag == -1
            alignment = 'right';
            NormShockTxtXL = NormShockXL;
        elseif HFlag == 1
            alignment = 'left';
            NormShockTxtXL=(xSL + 0.25*(xSR - xSL)+abs(min(xlim)))/diff(xlim)*pos(3) + pos(1);
        end
        
        SharrowX1X2 = [NormShockXL,NormShockXR];
        SharrowY1Y2 = [NormShockY,NormShockY];
        ShockTxtBox = [NormShockTxtXL, NormShockYTxt ,NormShockYTxtWdth, NormShockYTxtHgt];
        DeltaX = diff(SharrowX1X2);
        if DeltaX < 0.05
            ArrowSz = floor(DeltaX/0.05*10) + 1;
        else
            ArrowSz = 10;
        end
        annoShock = annotation('doublearrow',SharrowX1X2,SharrowY1Y2,'Linewidth',2,'Head1length',ArrowSz,'Head1Width', ArrowSz,'Head2length',ArrowSz,'Head2Width', ArrowSz);
        ProfileSplit = strsplit(Profile,'_');
        ProfileNameNonLaTex = ProfileSplit{1};
%         for stridx = 2:( length(ProfileSplit) )
%             ProfileNameNonLaTex = [ProfileNameNonLaTex,'-',ProfileSplit{stridx}];
%         end
        
        %Place 'profile name  width = length' along arrow that marks width of profile
        ShockTxt = sprintf('%s width = %4.1f', ProfileNameNonLaTex,LShock);
        
       
        annotation('textbox',ShockTxtBox,'EdgeColor','none','String',ShockTxt,'FontSize',13,'Linewidth',2,'HorizontalAlignment',alignment,'FitBoxToText','on');

        
        
        %*********************************************
        %% Plot vertical lines to Left and Right of Shock 
        % & horizontal line connecting them
        %*********************************************
%         plot([xSL xSR],[0.8*Zmin, 0.8*Zmin],'Color',[0 0 0]); % Horizontal line underneath Shock
%         hold('on');
        
        plot([xSL xSL],[0.5*Zmin, 0.8*Zmax],'Color',[0 0 0]); % Vertical Line on LHS of Shock
        hold('on');
        
        plot([xSR xSR],[0.5*Zmin, 0.8*Zmax],'Color',[0 0 0]); % Vertical Line on RHS of Shock
        hold('on');

        %*********************************************
        %% Mark LHS of Shock
        %*********************************************
        
        XBzL = x(idSL); YBzL = imag(b(idSL));             % 'X' marker intersection of LHS with Im(b)
        plot(XBzL,YBzL,'x','MarkerFaceColor','red','MarkerSize',5);
        hold('on');

        idBzTxtL = idSL - floor(0.2*idOS);              % Text label, value of Im(b) at LHS
        if idBzTxtL < 1
            idBzTxtL = 1;
            TxtLAlign = 'left';
        else
            TxtLAlign = 'right';
        end
        XBzTxtL = x(idBzTxtL); YBzTxtL = imag(b(idSL))-0.1*Zmax;
        text(XBzTxtL,YBzTxtL,sprintf('Im(b) = %s ',num2str(imag(b(idSL)))),'Color',[1 0 0],'horizontalAlignment', TxtLAlign);
        hold 'on';

        XBmagL = x(idSL); YBmagL = abs(b(idSL));          % 'X' marker intersection of LHS with Abs(b)
        plot(XBmagL,YBmagL,'x','MarkerFaceColor','blue','MarkerSize',5);
        hold 'on';

        idBmagTxtL = idSL - floor(0.2*idOS);            % Text label, value of Abs(b) at LHS
        if idBmagTxtL < 1
            idBmagTxtL = 1;
            TxtLAlign = 'left';
        else
            TxtLAlign = 'right';
        end
        XBmagTxtL = x(idSL); YBzTxtL = abs(b(idSL))+ 0.1*Zmax;
        text(XBmagTxtL,YBzTxtL,sprintf('Abs(b) = %s ',num2str(abs(b(idSL)))),'Color',[0 0 1],'horizontalAlignment', TxtLAlign);
        hold 'on';

        %*********************************************
        %% Mark RHS of Shock
        %*********************************************

        if contains(Profile, 'HybridIS')
            YTxtBaseRHS = Zmax;
        else
            YTxtBaseRHS = imag(b(idSR));
        end
        
        XBzR = x(idSR); YBzR = imag(b(idSR));
        plot(XBzR,YBzR,'x','MarkerFaceColor','red','MarkerSize',5); % 'X' marker intersection of RHS with Im(b)
        hold 'on';

        idBzTxtR = idSR + floor(0.2*idOS);              % Text label, value of Im(b) at RHS
        if idBzTxtR > Nmax
            idBzTxtR = Nmax;
            TxtRAlign = 'right';
        else
            TxtRAlign = 'left';
        end
        XBzTxtR = x(idBzTxtR); YBzTxtR = imag(b(idSR))-0.1*Zmax;
        text(XBzTxtR,YBzTxtR,sprintf('Im(b) = %s ',num2str(imag(b(idSR)))),'Color',[1 0 0],'horizontalAlignment', TxtRAlign);
        hold 'on';

        XBmagR = x(idSR); YBmagR = abs(b(idSR));          % 'X' marker intersection of RHS with Abs(b)
        plot(XBmagR,YBmagR,'x','MarkerFaceColor','blue','MarkerSize',5);
        hold 'on';

        idBmagTxtR = idSR + floor(0.2*idOS);            % Text label, value of Abs(b) at RHS
        if idBmagTxtR > Nmax
            idBmagTxtR = Nmax;
            TxtRAlign = 'right';
        else
            TxtRAlign = 'left';
        end
        XBmagTxtR = x(idBmagTxtR); YBmagTxtR = abs(b(idSL))+ 0.1*Zmax;
        text(XBmagTxtR,YBmagTxtR,sprintf('Abs(b) = %s ',num2str(abs(b(idSR)))),'Color',[0 0 1],'horizontalAlignment', TxtRAlign);
        hold 'on';
        
        idCRTxtR = idSR + floor(0.2*idOS);            % Text label, value of Abs(b) at RHS
        if idCRTxtR > Nmax
            idCRTxtR = Nmax;
            TxtRAlign = 'right';
        else
            TxtRAlign = 'left';
        end
        XCRTxtR = x(idCRTxtR); YCRTxtR = abs(b(idSL)) - 0.1*Zmax;
        text(XCRTxtR,YCRTxtR,sprintf('|Compression Ratio| = %s ',num2str( abs(b(idSR))/abs(b(idSL)) )),'Color',[0 0 1],'horizontalAlignment', TxtRAlign);
        hold 'on';

    end
    
    function sname = getname(argName)
        sname = inputname(argName);
    end
    function [Bflag, FP, Prfl, xLprofile, xRprofile, Lprofile, idLnew, idRnew] = Boundary_Approval(Fig, Profile, Prfl, prflKey, FP, Xp, XSkL, XSkR)
        figure(Fig);
        xLprofile = XSkL;
        xRprofile = XSkR;
        Lprofile = (XSkR - XSkL);
        Prfl(prflKey) = Lprofile;
        idLnew = 1; idRnew = length(Xp);
        Bflag = 1;
        if FP == 1 
            FP = FP + 1;
            % Construct a questdlg with three options
            choiceTxt = sprintf('Accept suggested %s boundaries?',Profile);
            choice = questdlg(choiceTxt, ...
                'Accept marked boundaries?', ...
                'Yes','No','No');
            % Handle response
            switch choice
                case 'Yes'
                    disp([choice ' Woot!'])
                    Lprofile = (XSkR - XSkL);
                    Prfl(prflKey) = Lprofile;
%                     idLsol = idxShockL;
%                     idRsol = idxShockR;
                    xLprofile = XSkL;
                    xRprofile = XSkR;
                    tmpL = abs(Xp - xLprofile);
                    tmpR = abs(Xp - xRprofile);
                    [val, idLnew] = min(tmpL); [val, idRnew] = min(tmpR); 
                    idLnew = round(idLnew/2)*2;
                    Bflag = 2;
                case 'No'
                    disp([choice ' Overriding'])
                    Bflag = 1;
                case 'No'
                    disp([choice ' Overriding'])
                    Bflag = 1;
            end
        else
            [x,y] = ginput(2);
            x = sort(x);
            xL = x(1); xR = x(2);
            tmpL = abs(Xp - xL); tmpR = abs(Xp - xR);
            [val, idLnew] = min(tmpL); [val, idRnew] = min(tmpR); 
            idLnew = round(idLnew/2)*2;

            % Construct a questdlg with three options
            choice = questdlg('Selection correct?', ...
                'Select Left-Right boundary of soliton:', ...
                'Yes','Zoom Out','Zoom In','Get Me Outa Here');
            % Handle response
            switch choice
                case 'Zoom Out'
                    disp([choice '! Zooming Out Now.'])
                    xLprofile = Xp(1);
                    xRprofile = Xp(end);
                    Lprofile = xRprofile - xLprofile;
                    idLnew = 1;
                    idRnew = length(Xp);
                    Bflag = 1;
                case 'Yes'
                    disp([choice ' Woot!'])
                    tmpL1 = abs(Xp - xL);
                    tmpR1 = abs(Xp - xR);
                    [val, idLnew] = min(tmpL1);
                    [val, idRnew] = min(tmpR1);
                    Lprofile = (xR - xL);
                    Prfl(prflKey) = Lprofile;
%                     idLsol = idL1;
%                     idRsol = idR1;
                    xLprofile = xL;
                    xRprofile = xR;
                    Bflag = 2;
                case 'Zoom In'
                    disp([choice ' Zoom!'])
                    tmpL1 = abs(Xp - xL);
                    tmpR1 = abs(Xp - xR);

                    [val, idLnew] = min(tmpL1);
                    [val, idRnew] = min(tmpR1);

                    Bflag = 1;
                case 'Get Me Outa Here'
                    disp('Prepare for Transport.')
                    Bflag = -1;
            end
        end
    end