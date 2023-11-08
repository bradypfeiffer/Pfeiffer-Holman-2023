function Real_EV_Excel_Archive(PStr, ReEV_List, time)
%% Archiving parameters stored in structure, PStr
% Called from Call_DNLS_v8.m right after call to DNLS_v4.m
% PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
%     'GCon',Graphics,'FCon',Folder,'DCon',Diagnostics);
    TimeExcelArchiveStart = cputime;
    FCon = PStr.FCon;
    TCon = PStr.TCon;
    NCon = PStr.NCon;
    dt   = NCon('dt');
    timestepsPerFrame = NCon('timestepsPerFrame');
    EVFolderAbs = [FCon('ArchiveBase'),'\Trial ',TCon('TrialNumber'),'\EV'];
    utilities(EVFolderAbs,1);
    SheetName = ['Lines_',num2str(ReEV_List{1}(1)),'_to_',num2str(ReEV_List{end}(1))];
    xlsfile = [EVFolderAbs,'\ReEVs_Trial_',TCon('TrialNumber'),'_Run_',TCon('RunNumber'),'.xlsx'];
    try
        %Check if an Excel server is running
        ex = actxGetRunningServer('Excel.Application');
    catch ME
        disp(ME.message)
    end
    if exist('ex','var')
        % Get the names of all open Excel files
        wbs = ex.Workbooks;
        % List the entire path of all excel workbooks that are currently open
        % then close each of them ...
        for i = 1:wbs.Count
            wbs.Item(i).FullName
            ['Saving ',wbs.Item(i).FullName]
            wbs.Item(i).Saved = 1;
            ['Closing ',wbs.Item(i).FullName]
            wbs.Item(i).Close;
            
        end
    end

    % https://www.mathworks.com/matlabcentral/answers/351328-how-to-save-and-close-excel-file-using-actxserver
    % https://www.mathworks.com/matlabcentral/answers/91547-how-do-i-insert-a-matlab-figure-into-an-excel-workbook-through-activex?s_tid=srchtitle
    e = actxserver('Excel.Application');
    File=xlsfile;
    Sheet = SheetName;
    % % Create Excel workbook if it does not already exist
    if ~exist(File,'file')
        excelWb = e.workbooks.Add;
        excelWb.SaveAs(File);
        excelWb.Close();
    end

    % Open Excel workbook
    invoke(e.Workbooks,'Open',File);

    Trial = PStr.TCon;  % used to print updates to Command Window
    
    %% Write Profile Parameters to archive
    Lcell = 'A'; Rcell = 'B'; topPrfl = 3;
    Params = PStr.PrflCon;
    heading = Params('title');
    elapsed = num2str(cputime - time); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell, Rcell, topPrfl,Params, heading, e)
    
    %% Write Graphics Parameters to archive
    Lcell = 'A'; Rcell = 'B'; topG = topPrfl + length(Params.keys) + 2;
    Params = PStr.GCon;
    heading = Params('title');
    elapsed = num2str(cputime - time); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell, Rcell, topG,Params, heading, e)
         
    %% Write Physics Parameters
    Lcell = 'A'; Rcell = 'B'; topP = topG + length(Params.keys) + 2;
    Params = PStr.PCon;
    heading = Params('title');
    elapsed = num2str(cputime - t); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell,Rcell, topP,Params, heading,e)
      
    %% Write Numerics Parameters
    Lcell = 'A'; Rcell = 'B'; topN = topP + length(Params.keys) + 2;
    Params = PStr.NCon;
    heading = Params('title');
    elapsed = num2str(cputime - t); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell,Rcell, topN,Params, heading,e)
   
    %% Write Trial Parameters
    Lcell = 'D'; Rcell = 'E'; topN = 3;
    Params = PStr.TCon;
    heading = Params('title');
    elapsed = num2str(cputime - t); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell,Rcell, topN,Params, heading,e)
    
    %% Write Diagnostics
    Lcell = 'D'; Rcell = 'E'; topD = topN + length(Params.keys) + 2;
    Params = PStr.DCon;
    heading = Params('title');
    elapsed = num2str(cputime - t); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell,Rcell, topD,Params, heading,e)

    %% Write Folder Parameters
    Lcell = 'D'; Rcell = 'E'; topF = topD + length(Params.keys) + 2;
    Params = PStr.FCon;
    heading = Params('title');
    elapsed = num2str(cputime - t); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile, Sheet, Lcell, Rcell, topF, Params, heading, e)

    %% Write Real Eigenvalues

    % Determine NumberOfLines and NumberEvs = # of real eigenvalues + 1 (for line #)
    NumberEvs = 1;
    NumberOfLines = length(ReEV_List);
    for EVidx = 1:NumberOfLines
        EV_length_current = length(ReEV_List{EVidx});
        if NumberEvs < EV_length_current
            NumberEvs = EV_length_current;
        end        
    end
   
    EV_Table = table;
    EV_entry = cell(NumberOfLines,NumberEvs);
    EV_entry(:,:)={0};
    for EVidx = 1:NumberOfLines
        
        EVs_current = ReEV_List{EVidx};
        for EVNidx = 1:length(EVs_current)   
            EV_entry(EVidx,EVNidx)= {ReEV_List{EVidx}(EVNidx)};
        end
    end

    Rcell = ['G','H','I','J','K','L','M','N','O','P'];    

    EV_Table_VariableNames{1} = 'Line Number';
    EV_Table_VariableNames{2} = 'Energy';
    EV_Table_VariableNames{3} = 'FELandau';
    
    for EVidx = 4:NumberEvs
        lambdan = ['lambda',num2str(EVidx - 3)];
        EV_Table_VariableNames{EVidx}=lambdan;
    end

    NumberEvs = NumberEvs + 1; % Add time as column (end + 1)
    EV_Table_VariableNames{NumberEvs} = 'Time';

    % Write Column Header Labels
    for EVidx=1:NumberEvs
        xlswrite1(xlsfile,EV_Table_VariableNames(EVidx), Sheet,[Rcell(EVidx),num2str(topD+0)],e);
    end

    % Write Real Eigenvalues by line
    for LNidx = 1:NumberOfLines
        TimeLN = (EV_entry{LNidx,1} - 1)*timestepsPerFrame*dt;
        xlswrite1(xlsfile,EV_entry{LNidx,1}, Sheet,[Rcell(1),num2str(topD+LNidx)],e);
       
        for EVidx=1:(NumberEvs - 1)
            xlswrite1(xlsfile,EV_entry{LNidx,EVidx}, Sheet,[Rcell(EVidx),num2str(topD+LNidx)],e);
        end
         xlswrite1(xlsfile,{TimeLN}, Sheet,[Rcell(NumberEvs),num2str(topD+LNidx)],e);
    end

    %BOLD Column Header
    EVheaderRange = [Rcell(1),num2str(topD),':',Rcell(NumberEvs),num2str(topD)];
    ec = e.ActiveSheet.Range(EVheaderRange);
    set(ec.Font, 'Bold',true)
    
    %Border Around Header and Data
    whole_range = [Rcell(1),num2str(topD),':',Rcell(NumberEvs),num2str(topD + NumberOfLines)];
    rngObj = e.ActiveSheet.Range(whole_range);
    set(rngObj.Borders,'LineStyle',9);

    % ExAct = e.Activesheet;   
    % ExActRange = get(ExAct,'Range',range);
    % ExActRange.Select;
    invoke(rngObj.Columns,'Autofit');

    


    ['Excel Archiving Completed: ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    ['Total Excel Archive Time = ',num2str(cputime - TimeExcelArchiveStart), ' s']
    invoke(e.ActiveWorkbook,'Save');

    % Close Excel workbook and delete and clear instance
    e.Quit;
    e.delete;
    clear e;  


end