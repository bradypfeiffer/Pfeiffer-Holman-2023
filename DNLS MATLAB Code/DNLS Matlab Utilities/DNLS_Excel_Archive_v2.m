function DNLS_Excel_Archive_v2(xlsfile, PStr, k, time)
%% Archiving parameters stored in structure, PStr
% Called from Call_DNLS_v8.m right after call to DNLS_v4.m
% PStr = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
%     'GCon',Graphics,'FCon',Folder,'DCon',Diagnostics);
    TimeExcelArchiveStart = cputime;
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


    e = actxserver ('excel.application');
    File=xlsfile;
    Sheet = strcat('Sheet',int2str(k));
    % Create Excel workbook if it does not already exist
    if ~exist(File,'file')
        excelWb = e.workbooks.Add;
        excelWb.SaveAs(File);
        excelWb.Close(false);
    end
    % Open Excel workbook
    invoke(e.Workbooks,'Open',File);

    Trial = PStr.TCon;  % used to print updates to Command Window
    
    %Write Profile Parameters to archive
    Lcell = 'A'; Rcell = 'B'; topPrfl = 3;
    Params = PStr.PrflCon;
    heading = Params('title');
    elapsed = num2str(cputime - time); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell, Rcell, topPrfl,Params, heading, e)
    
    %Write Graphics Parameters to archive
    Lcell = 'A'; Rcell = 'B'; topG = topPrfl + length(Params.keys) + 2;
    Params = PStr.GCon;
    heading = Params('title');
    elapsed = num2str(cputime - time); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell, Rcell, topG,Params, heading, e)
         
    %Write Physics Parameters
    Lcell = 'A'; Rcell = 'B'; topP = topG + length(Params.keys) + 2;
    Params = PStr.PCon;
    heading = Params('title');
    elapsed = num2str(cputime - t); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell,Rcell, topP,Params, heading,e)
      
    %Write Numerics Parameters
    Lcell = 'A'; Rcell = 'B'; topN = topP + length(Params.keys) + 2;
    Params = PStr.NCon;
    heading = Params('title');
    elapsed = num2str(cputime - t); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell,Rcell, topN,Params, heading,e)
   
    % Write Trial Parameters
    Lcell = 'D'; Rcell = 'E'; topN = 3;
    Params = PStr.TCon;
    heading = Params('title');
    elapsed = num2str(cputime - t); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell,Rcell, topN,Params, heading,e)
    
    %Write Diagnostics
    Lcell = 'D'; Rcell = 'E'; topD = topN + length(Params.keys) + 2;
    Params = PStr.DCon;
    heading = Params('title');
    elapsed = num2str(cputime - t); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell,Rcell, topD,Params, heading,e)

    %Write Folder Parameters
    Lcell = 'D'; Rcell = 'E'; topF = topD + length(Params.keys) + 2;
    Params = PStr.FCon;
    heading = Params('title');
    elapsed = num2str(cputime - t); t = cputime;
    ['Writing ',heading, ' to Excel for ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    Arrayxlswrite(xlsfile,Sheet,Lcell,Rcell, topF,Params, heading,e)
    ['Excel Archiving Completed: ', Trial('TrialName'), ' --- Elapsed time = ',elapsed, ' s']
    ['Total Excel Archive Time = ',num2str(cputime - TimeExcelArchiveStart), ' s']
    invoke(e.ActiveWorkbook,'Save');

    % Close Excel workbook and delete and clear instance
    e.Quit;
    e.delete;
    clear e;  


end