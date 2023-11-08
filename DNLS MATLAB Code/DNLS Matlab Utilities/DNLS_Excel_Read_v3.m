function [PStr] = DNLS_Excel_Read_v3(XLfile, sheet)
%% Comments
%Modified 6/12/2016
%If additional Key/Value pairs are added later they can be read,
%but ONLY IF the 'Key' variables are passed to DNLS_Excel_Read
%The advantage of this version of the code is that having additional Key
%variables and values will not cause a run time error or, worse, improperly
%evaluate a Key variable. 
%% _Function Definitions
function DumValues = Values(DumKey, DumTitle, profile)
    DTitleDebug = 'nada';
    if strcmp(DumTitle, DTitleDebug)
        ['*****************************************']
        ['     Starting Numerics Values Degugging:']
    end
    for ii = 1:length(DumKey)
        key = DumKey(ii);
        
        if strcmp(DumTitle, DTitleDebug)
            fprintf('top of loop');
            ii
            key
        end
                     
        if ~strcmp(key,'title') && ~strcmp(key,'Profile')
            RN = strcmp(raw(:,1),key); %returns row #(s) in which 'Key' appears
            temp1 = raw(find(RN ==1),2); %Assuming that Key appears just once
            if strcmp(DumTitle, DTitleDebug)
                fprintf('if key NOT title or Profile');i
                fprintf('temp1 = ');
                temp1
                ['RN = ']
                RN
                ii
                key
            end
        elseif strcmp(key,'title')
            temp1 = {DumTitle};
            if strcmp(DumTitle, DTitleDebug)
                fprintf('if key IS title');
                ii
                key
            end
        elseif strcmp(key,'Profile')
            temp1 = {profile};           
            if strcmp(DumTitle, DTitleDebug)
                fprintf('if key IS Profile');
                ii
                key
            end
        end
        if isempty(temp1)
            temp1 = {'NA'};
            if strcmp(DumTitle, DTitleDebug)
                fprintf('temp1 is empty so set temp1 = %s',temp1);
                ii
                key
            end
        end
        DumValues(ii) = temp1(1);
    end
end
%% Read Excel Archive            
Sheet = ['Sheet',num2str(sheet)];
 %Post Run Diagnostics
    [num,txt,rawDE] = xlsread(XLfile,Sheet,'D3:E60');
    [num,txt,rawAB] = xlsread(XLfile,Sheet,'A3:B60');
    raw = [rawAB;rawDE];
%% Determine profile type    
%Read profile type. If not found set as '1psoliton'
    RN = strcmp(raw(:,1),'Profile'); %returns row # in which 'Profile' appears
    Profile = '1psoliton'; %default: earliest profile, prior to there being two or more options
    if max(RN) == 1
        ProfileTemp = raw(find(RN ==1),2); %possibly Helical, 1psoliton, or newer profile
        Profile = ProfileTemp{:};
    end
%% Determine Numerical Method used    
%Read profile type. If not found set as '1psoliton'
    RNTN = strcmp(raw(:,1),'TrialName'); %returns row # in which 'Profile' appears
    Profile = 'DNLS_v1'; %default: earliest profile, prior to there being two or more options
    if max(RNTN) == 1
        ProfileTemp = raw(find(RN ==1),2); %possibly Helical, 1psoliton, or newer profile
        Profile = ProfileTemp{:};
    end
%% Set Prfl Keys based on profile type    
    %Set PhysicsKeys based on profile type
    %update Prfl Keys below as new profiles are implemented
    PrflKeys   = {'title', 'Profile', 'speed'}; %default
    if strcmp(Profile,'Helical')
        PrflKeys = {'title', 'Profile','speed', 'W_AS', 'W_0', 'L_AS', 'K_Helix', 'helix_angle', 'theta_p', 'xL', 'xR','EVL', 'EVR'};
    end  
    if strcmp(Profile,'1psoliton')
        PrflKeys   = {'title', 'Profile', 'speed', 'lambda',  'sign', 'bl'};
    end
    if strcmp(Profile,'2psoliton')
        PrflKeys   = {'title', 'Profile', 'speed', 'lr',  'li'};
    end
    if strcmp(Profile,'FastShock')
        PrflKeys   = {'title'   , 'Profile', 'speed', 'bl',  'br', 'xfs', 'xrhs', 'wfs', 'wrhs'};
    end
%% Define Keys for containers    
    Ttitle = 'Trial Identifiers';
    TrialKeys = {'title', 'RunNumber','TrialNumber','TrialName','TimeStamp','TrialObjective','TrialID'};
    TrialValues = cell(1,length(TrialKeys));
    
    Ntitle = 'Numerics Parameters';
    NumericsKeys = {'title', 'N', 'L', 'dx', 'Courant', 'dt', 'timestepsPerFrame', 'TimestepsAttempted'};
    NumericsValues = cell(1,length(NumericsKeys));
   
    Ptitle = 'Physics Parameters';
    PhysicsKeys   = {'title','bo','nu'};
    PhysicsValues = cell(1,length(PhysicsKeys)); 
    
    Prfltitle = 'Profile Parameters';
    PrflValues = cell(1,length(PrflKeys));
    
    Gtitle = 'Graphics Parameters';
    GraphicsKeys = {'title', 'DisplayExactSolution','zmax','zmin'};
    GraphicsValues = cell(1,length(GraphicsKeys));
    
    Ftitle = 'File Locations';
    FolderKeys = {'title', 'ArchiveBase', 'TrialFolder', 'BfieldFolder','Bfieldfile'};
    FolderValues = cell(1,length(FolderKeys)); 
    
    Dtitle = 'Post Run Diagnostics';
    DiagnosticKeys = {'title', 'CpuTimeUsed', 'timesteps','growth', 'decay'};
    DiagnosticValues = cell(1,length(DiagnosticKeys)); 
%%Search raw excel data for values   
TrialValues      = Values(TrialKeys, Ttitle, Profile);
NumericsValues   = Values(NumericsKeys, Ntitle, Profile);
% ['From Excel_Read_v3: <Search raw excel data for values>']
% NumericsValues
PhysicsValues    = Values(PhysicsKeys, Ptitle, Profile);
PrflValues    = Values(PrflKeys, Prfltitle, Profile);
GraphicsValues   = Values(GraphicsKeys, Gtitle, Profile);
FolderValues     = Values(FolderKeys, Ftitle, Profile);
DiagnosticValues = Values(DiagnosticKeys, Dtitle, Profile);

%%Prepare Containers
Trial = containers.Map(TrialKeys, TrialValues);
Numerics = containers.Map(NumericsKeys, NumericsValues);
% ['From Excel_Read_v3: <Prepare Containers>']
% NumericsValues
% Numerics.values
Physics = containers.Map(PhysicsKeys, PhysicsValues);
Prfl = containers.Map(PrflKeys, PrflValues);
Graphics = containers.Map(GraphicsKeys, GraphicsValues);
Folder = containers.Map(FolderKeys, FolderValues);
Diagnostics = containers.Map(DiagnosticKeys, DiagnosticValues);

%%Finish Parameter Structure with categorized containers (keys v. values)
PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
    'PrflCon',Prfl,'GCon',Graphics,'FCon',Folder,'DCon',Diagnostics);
end