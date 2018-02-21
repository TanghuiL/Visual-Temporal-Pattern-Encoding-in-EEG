% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % Generating trigger codes from data (xml) exported from EXPO% % % %
% % % % % % % % % % VisA study (for Active Condition)% % %  % % %  % %
% % % % % % % % % % % % Dr Sussman's lab % % % % % % % % % % % % % % %
% % % % % % % % % Albert Einstein College of Medicine % % % % % % % % %
% % % % Last updated on 07/27/2016 by Huizhen Tang (Joann) % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Get xml files
clear
[filename,pathname] = uigetfile(...
    { '*.xml','EXPO exported data files (*.xml)';'*.*','All Files' }, ...
    'Select xml file(s)','Multiselect','on');
% Abort if the user hit 'Cancel'
if isequal(filename,0)||isequal(pathname,0),
    disp('Aborted.');
    return;
end
%% This line is very important! The xml files need to be in the current folder!
cd(pathname)
%% getting indicator for differentiating trigger for different conditions 
prompt = {'Please indicate the condition being processed (1 - 10deg; 2 - 90deg)'};
dlg_title = 'Which condition?';
num_lines = 1;
defAns = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,defAns);%If the user clicks the Cancel button to close an input dialog box,
% % % Abort if the user clicks 'Cancel'.
if isempty(answer), disp('Aborted.');
    return;
end
[cond status] = str2num(answer{1});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid ','Error in Parameter settings','error');
end
%% extracting codes from EXPO xml file % % % % % % % % % % % % %
attriALL = [];
if iscell(filename)
    for i = 1:length(filename)
        [s] = xml2struct(filename{i});%converting xml data to struct
        % Reading the attributs 
        m = 0; block = []; slot = []; strTime = []; endTime = [];
        for n = 1:size(s.ExpoXData.Passes.Pass,2)
            BlockID = s.ExpoXData.Passes.Pass{1,n}.Attributes.BlockID;
            bID = str2num(BlockID);
            if bID == 0 || bID == 3,
                m = m +1;
            else
                k = n - m;
                block(k) = bID;
                SlotID = s.ExpoXData.Passes.Pass{1,n}.Attributes.SlotID;
                slot(k) = str2num(SlotID);
                StartTime = s.ExpoXData.Passes.Pass{1,n}.Attributes.StartTime;
                strTime(k) = str2num(StartTime);
                EndTime = s.ExpoXData.Passes.Pass{1,n}.Attributes.EndTime;
                endTime(k) = str2num(EndTime);
            end
        end
        attri = cat(1,block,slot,strTime,endTime);
        attriALL = attri';
       % recoding trigger using slotID to differentiate stimulus's codition as well as its position
        trgg = []; mm = 0;
        if cond == 1 %%% condition "10degree"
            for nn = 1:size(attriALL,1)
                slot = attriALL(nn,2);
                % % % % % % recoding for Deviant 1 pattern (XXO) % % % % % % % %
                if slot == 16
                    kk = nn - mm; trgg(kk) = 41;
                elseif slot == 18
                    kk = nn - mm; trgg(kk) = 42;
                elseif slot == 20
                    kk = nn - mm; trgg(kk) = 40;
                    % % % % % % recoding for Deviant 2 pattern (XXXXO) % % % % % % % %
                elseif slot == 22
                    kk = nn - mm; trgg(kk) = 51;
                elseif slot == 24
                    kk = nn - mm; trgg(kk) = 52;
                elseif slot == 26
                    kk = nn - mm; trgg(kk) = 53;
                elseif slot == 28
                    kk = nn - mm; trgg(kk) = 54;
                elseif slot == 30
                    kk = nn - mm; trgg(kk) = 50;
                    % % % % % % recoding for Standard pattern (XXXO) % % % % % % % %
                elseif slot == 8
                    kk = nn - mm; trgg(kk) = 11;
                elseif slot == 10
                    kk = nn - mm; trgg(kk) = 12;
                elseif slot == 12
                    kk = nn - mm; trgg(kk) = 13;
                elseif slot == 14
                    kk = nn - mm; trgg(kk) = 10;
                else
                    mm = mm +1;
                end
            end
        elseif cond == 2 %%% condition "90degree"
            for nn = 1:size(attriALL,1)
                slot = attriALL(nn,2);
                % % % % % % recoding for Deviant 1 pattern (XXO) % % % % % % % %
                if slot == 16
                    kk = nn - mm; trgg(kk) = 141;
                elseif slot == 18
                    kk = nn - mm; trgg(kk) = 142;
                elseif slot == 20
                    kk = nn - mm; trgg(kk) = 140;
                    % % % % % % recoding for Deviant 2 pattern (XXXXO) % % % % % % % %
                elseif slot == 22
                    kk = nn - mm; trgg(kk) = 151;
                elseif slot == 24
                    kk = nn - mm; trgg(kk) = 152;
                elseif slot == 26
                    kk = nn - mm; trgg(kk) = 153;
                elseif slot == 28
                    kk = nn - mm; trgg(kk) = 154;
                elseif slot == 30
                    kk = nn - mm; trgg(kk) = 150;
                    % % % % % % recoding for Standard pattern (XXXO) % % % % % % % %
                elseif slot == 8
                    kk = nn - mm; trgg(kk) = 111;
                elseif slot == 10
                    kk = nn - mm; trgg(kk) = 112;
                elseif slot == 12
                    kk = nn - mm; trgg(kk) = 113;
                elseif slot == 14
                    kk = nn - mm; trgg(kk) = 110;
                else
                    mm = mm +1;
                end
            end
        else
            disp('Condition does not exist')
        end
        trgg = trgg';
        save([pathname filename{i}(1:(end-4))],'attriALL','trgg')
    end
else
    [s] = xml2struct(filename);%converting xml data to struct
        %%%% Reading the attributs and put values in different run to one variable
        m = 0; block = []; slot = []; strTime = []; endTime = []; 
        for n = 1:size(s.ExpoXData.Passes.Pass,2)
            BlockID = s.ExpoXData.Passes.Pass{1,n}.Attributes.BlockID;
            bID = str2num(BlockID);
            if bID == 0 || bID == 3,
                m = m +1;
            else
                k = n - m;
                block(k) = bID;
                SlotID = s.ExpoXData.Passes.Pass{1,n}.Attributes.SlotID;
                slot(k) = str2num(SlotID);
                StartTime = s.ExpoXData.Passes.Pass{1,n}.Attributes.StartTime;
                strTime(k) = str2num(StartTime);
                EndTime = s.ExpoXData.Passes.Pass{1,n}.Attributes.EndTime;
                endTime(k) = str2num(EndTime);
                
            end
        end
        attri = cat(1,block,slot,strTime,endTime);
        attriALL = attri';
        %%%% recoding trigger using slotID to differentiate stimulus's codition as well as its position
        trgg = []; mm = 0;
        if cond == 1 %%% condition "10degree"
            for nn = 1:size(attriALL,1)
                slot = attriALL(nn,2);
                % % % % % % recoding for Deviant 1 pattern (XXO) % % % % % % % %
                if slot == 16
                    kk = nn - mm; trgg(kk) = 41;
                elseif slot == 18
                    kk = nn - mm; trgg(kk) = 42;
                elseif slot == 20
                    kk = nn - mm; trgg(kk) = 40;
                    % % % % % % recoding for Deviant 2 pattern (XXXXO) % % % % % % % %
                elseif slot == 22
                    kk = nn - mm; trgg(kk) = 51;
                elseif slot == 24
                    kk = nn - mm; trgg(kk) = 52;
                elseif slot == 26
                    kk = nn - mm; trgg(kk) = 53;
                elseif slot == 28
                    kk = nn - mm; trgg(kk) = 54;
                elseif slot == 30
                    kk = nn - mm; trgg(kk) = 50;
                    % % % % % % recoding for Standard pattern (XXXO) % % % % % % % %
                elseif slot == 8
                    kk = nn - mm; trgg(kk) = 11;
                elseif slot == 10
                    kk = nn - mm; trgg(kk) = 12;
                elseif slot == 12
                    kk = nn - mm; trgg(kk) = 13;
                elseif slot == 14
                    kk = nn - mm; trgg(kk) = 10;
                else
                    mm = mm +1;
                end
            end
        elseif cond == 2 %%% condition "90degree"
            for nn = 1:size(attriALL,1)
                slot = attriALL(nn,2);
                % % % % % % recoding for Deviant 1 pattern (XXO) % % % % % % % %
                if slot == 16
                    kk = nn - mm; trgg(kk) = 141;
                elseif slot == 18
                    kk = nn - mm; trgg(kk) = 142;
                elseif slot == 20
                    kk = nn - mm; trgg(kk) = 140;
                    % % % % % % recoding for Deviant 2 pattern (XXXXO) % % % % % % % %
                elseif slot == 22
                    kk = nn - mm; trgg(kk) = 151;
                elseif slot == 24
                    kk = nn - mm; trgg(kk) = 152;
                elseif slot == 26
                    kk = nn - mm; trgg(kk) = 153;
                elseif slot == 28
                    kk = nn - mm; trgg(kk) = 154;
                elseif slot == 30
                    kk = nn - mm; trgg(kk) = 150;
                    % % % % % % recoding for Standard pattern (XXXO) % % % % % % % %
                elseif slot == 8
                    kk = nn - mm; trgg(kk) = 111;
                elseif slot == 10
                    kk = nn - mm; trgg(kk) = 112;
                elseif slot == 12
                    kk = nn - mm; trgg(kk) = 113;
                elseif slot == 14
                    kk = nn - mm; trgg(kk) = 110;
                else
                    mm = mm +1;
                end
            end
        else
            disp('Condition does not exist')
        end
        trgg = trgg';
        save([pathname filename(1:(end-4))],'attriALL','trgg')
end