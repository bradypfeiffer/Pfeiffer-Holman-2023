function utilities(FolderAbs, FolderRel)

if ~exist(FolderAbs, 'dir')
   % [FolderRel, ' not found. Will create it now.']
    %choice = menu([TrialFolder, ' not found.'],'I will make it now');
    if ~mkdir(FolderAbs)
        %choice = menu([folder, ' not found.'],'I will make it now');
    else %[FolderRel, ' created!']
    end
else %[FolderRel, ' found!']
end
end