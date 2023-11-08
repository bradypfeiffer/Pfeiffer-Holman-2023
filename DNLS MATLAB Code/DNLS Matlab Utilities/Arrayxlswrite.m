function Arrayxlswrite(filename, Sheet, L, R, top, Params, heading,e)

%     Excel = evalin('base','Excel');
    Excel = e;

    keys = Params.keys;
    values = Params.values;
    Nentries = length(keys) - 1;
    KeySz = length(keys);
    topcell = top;
    Lcol = L;
    Rcol = R;
    A1 = strcat(Lcol,int2str(topcell));
    B1 = strcat(Rcol,int2str(topcell));
    Bf = strcat(Rcol,int2str(topcell + Nentries));
    AtoB = strcat(Lcol,':',Rcol);
    A1toB1 = strcat(Lcol,int2str(topcell),':',Rcol,int2str(topcell));
    A1toBf = strcat(Lcol,int2str(topcell),':',Rcol,int2str(topcell+Nentries));
    B2toBf = strcat(Rcol,int2str(topcell+1),':',Rcol,int2str(topcell+Nentries));
    xlswrite1(filename,{heading},Sheet,A1,Excel);
    RowPrint = 0;
    Maxstring = -1;
    SpecialAlignment = 0;
%     eS = e.ActiveSheet;
%     CrntColWidth = eS.Range(AtoB).ColumnWidth;
    BTObj = 'F3';
    for k = 1:KeySz
        key = keys(k);
        valuec = values{k};
        valuep = values(k);
        if iscellstr(valuep) && ~strcmp(key,'title')
            temp = strfind(valuec,'\');
            if length(temp) ~= 0
                index = temp(end) ;
                valuec = valuec(index:end);
            end
        end        
        if~strcmp(key,'title')
            RowPrint = RowPrint + 1;
            An = strcat(Lcol,int2str(RowPrint + topcell));
            Bn = strcat(Rcol,int2str(RowPrint + topcell));
            xlswrite1(filename,key     ,Sheet,An,Excel)
            xlswrite1(filename,{valuec},Sheet,Bn,Excel)
            if strcmp(key,'TrialObjective')
                SpecialAlignment = 1;
                BTObj = strcat(Bn,':',Bn);
                Maxstring = length(valuec);
            end
            if length(valuec) > Maxstring
                Maxstring = length(valuec);
            end
        end
    end
    
%     e = actxserver('excel.application');
%      eW = e.Workbooks;
     %eF = eW.Open(filename);
      eS = e.ActiveSheet;

% e.ActiveSheet.Range(AtoB).ColumnWidth = 25;
if Maxstring < 25
    Maxstring = 25;
end
% CrntColWidth = eS.Range(AtoB).ColumnWidth;
% if Maxstring < CrntColWidth
%     Maxstring = CrntColWidth;
% end
eSRange = eS.Range(AtoB);
eSRange.ColumnWidth = Maxstring;

%BOLD Column Header
ec = e.ActiveSheet.Range(A1toB1);
set(ec.Font, 'Bold',true)

dat_range = A1toBf; 
rngObj = eS.Range(dat_range);
set(rngObj.Borders,'LineStyle',9);

% eF.Save
% eF.Close
% e.Quit


xlsalign(e, filename,Sheet,A1toB1,'MergeCells',1);
xlsalign(e, filename,Sheet,A1toB1,'Horizontal',3);
xlsalign(e, filename,Sheet,B2toBf,'Horizontal',3);

if SpecialAlignment
    xlsalign(e, filename,Sheet,BTObj,'Horizontal',2);
end