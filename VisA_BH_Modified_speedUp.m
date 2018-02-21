% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % Behavioral analysis % % % % % % % % % % % %
% % % % % % % % % % VisA study (for Active Condition)% % %  % % %  % %
% % % % % % % % % % % % Dr Sussman's lab % % % % % % % % % % % % % % %
% % % % % % % % % Albert Einstein College of Medicine % % % % % % % % %
% % % % Last updated on Mar/02/2017 by Huizhen Tang (Joann) % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [attriALL] = VisA_BH_Modified_speedUp
clear
%% Get xml files
[filename,pathname] = uigetfile(...
    { '*.xml','EXPO exported data files (*.xml)';'*.*','All Files' }, ...
    'Select xml file(s)','Multiselect','on');
% Abort if the user hit 'Cancel'
if isequal(filename,0)||isequal(pathname,0),
    disp('Aborted.');
    return;
end
%% getting indicator for differentiating trigger for different conditions
prompt = {'Please indicate the condition being processed (eg TR = Task Relevant)', ...
    'Inter-pattern ISI (eg 500 ms)','Inter-grating ISI (eg 200 ms)', ...
    'Lower end of time window for button response (eg. 200 ms)', ...
    'Upper end of time window for button response (eg. 900 ms)'};
dlg_title = 'Experimental parameters';
num_lines = 1;
defAns = {'TR','350','350','200','900'};
answer = inputdlg(prompt,dlg_title,num_lines,defAns);%If the user clicks the Cancel button to close an input dialog box,
% Abort if the user clicks 'Cancel'.
if isempty(answer), disp('Aborted.');
    return;
end
cond = answer{1};
[ip_ISI status] = str2num(answer{2});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid input','Error in Parameter settings','error');
end
[ig_ISI status] = str2num(answer{3});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid input','Error in Parameter settings','error');
end
[tw_low status] = str2num(answer{4});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid input','Error in Parameter settings','error');
end
[tw_up status] = str2num(answer{5});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid input','Error in Parameter settings','error');
end
%% extracting codes from EXPO xml file
attriALL = [];
if iscell(filename)
    for i = 1:length(filename)
        [s] = xml2struct([pathname filename{i}]);%converting xml data to struct
        % Reading the attributs
        [attri] = readAttributes(s);
        attriALL = cat(1,attriALL,attri);
        %         save([pathname 'BH_' cond '_' filename{i}(1:(end-4))],'attri')
        disp(['Finished  ' filename{i}])
    end
    attA = attriALL;
else
    [s] = xml2struct([pathname filename]);%converting xml data to struct
    [attri] = readAttributes(s);
    %      save([pathname 'BH_' cond '_' filename(1:(end-4))],'attri')
    disp(['Finished ' filename])
    attA = attri;
end
%% Calculate behavioral indices
% nHit - hit rate; rFAR - false alarm; rt - reaction time
% nHit_d1; rt_d1 (when target is O in deviant pattern 1 -- XXO)
% nHit_d2_x; rt_d2_x (when target is X4 in deviant pattern 2 -- XXXXO)
% nHit_d2_o; rt_d2_o (when target is X4 in deviant pattern 2 -- XXXXO)
nHit_d1 = 0; nHit_d2_x = 0; nHit_d2_o = 0; nFar = 0;
rt_d1 = []; rt_d2_x = [];rt_d2_o = []; nRsp = 0;
nd1 = 0; nd2 = 0;
for nr = 1:size(attA,1),
    if ~isnan(attA(nr,4)),
        nRsp = nRsp + 1;
    end
end

if ip_ISI == 500
    resp_tw =  5;% 1300 ms = 200 + 500 + 200 + 200 + 200
elseif ip_ISI == 350
    resp_tw = 4 % 1355 ms = 350 + 350 + 350 + 350
elseif ip_ISI == 200
    resp_tw = 6; % 1200 ms = 200 + 200 + 200 + 200 + 200 + 200
end



for i = 1:size(attA,1)
    slot = attA(i,2);
    if slot == 20, % target is deviant 1 (XXO)
        nd1 = nd1 + 1;
        for j = i:i+resp_tw
            if ~isnan(attA(j,4)),
                nHit_d1 = nHit_d1 + 1;
                if attA(j,4)> attA(j,5) && attA(j,4) < attA(j,6)
                    tem_rt = attA(j,4) - attA(i,5);
                    if tem_rt > tw_low && tem_rt < tw_up
                        rt_d1(nHit_d1) = tem_rt;
                    else
                    end
                else
                    disp([num2str(j) '--invalid time'])
                end
            end
        end
    elseif slot == 28, % target is deviant 2 XXXXO
        nd2 = nd2 + 1;
        for j = i:i+resp_tw
            if ~isnan(attA(j,4)),
                nHit_d2_x = nHit_d2_x + 1;
                if attA(j,4)> attA(j,5) && attA(j,4) < attA(j,6)
                    tem_rt = attA(j,4) - attA(i,5);
                    if tem_rt  > tw_low && tem_rt < tw_up
                        rt_d2_x(nHit_d2_x) = tem_rt;
                    else
                    end
                else
                    disp([num2str(j) '--invalid time'])
                end
            end
        end
    elseif slot == 30, % target is deviant 2 XXXXO
        for j = i:i+resp_tw
            if ~isnan(attA(j,4)),
                nHit_d2_o = nHit_d2_o + 1;
                if attA(j,4)> attA(j,5) && attA(j,4) < attA(j,6)
                    tem_rt = attA(j,4) - attA(i,5);
                    if tem_rt > tw_low && tem_rt < tw_up
                        rt_d2_o(nHit_d2_o) = tem_rt;
                    else
                    end
                else
                    disp([num2str(j) '--invalid time'])
                end
            end
        end
    end
end

%remove all the responses with a reaction time < 300 ms or > 900 ms; which
%sets the time window for correct response as 300 - 900 ms
if ~isempty(rt_d1)
    rt_d1(rt_d1<300) = []; rt_d1(rt_d1>900) =[];
end
if ~isempty(rt_d2_x)
    rt_d2_x(rt_d2_x<300) = []; rt_d2_x(rt_d2_x>900) = [];
end
if ~isempty(rt_d2_o)
    rt_d2_o(rt_d2_o<300) = []; rt_d2_o(rt_d2_o>900) = [];
end
mRt_d1 = mean(rt_d1);
mRt_d2_x = mean(rt_d2_x);
mRt_d2_o = mean(rt_d2_o);
rFar = (nRsp-length(rt_d1)-length(rt_d2_x))/nRsp;
rHit_d1 = length(rt_d1)/nd1;
rHit_d2_x = length(rt_d2_x)/nd2;
rHit_d2_o = length(rt_d2_o)/nd2;
rHit = (length(rt_d1) + length(rt_d2_x))/(nd1+nd2);
save([pathname 'BH_' cond],'nd1','nd2','nRsp','rHit','rHit_d1','rHit_d2_x',...
    'rHit_d2_o','mRt_d1','mRt_d2_x','mRt_d2_o','rFar','attA','rt_d1','rt_d2_x','rt_d2_o')
end
%% function for reading attributes from the XML files
function [attri] = readAttributes(s) % s is a struct
% Reading the attributs
m = 0; block = []; slot = []; strT = []; endT = [];
rid = [];rspt = []; pass = [];
for n = 1:size(s.ExpoXData.Passes.Pass,2)
    SlotID = s.ExpoXData.Passes.Pass{1,n}.Attributes.SlotID;
    PassID = s.ExpoXData.Passes.Pass{1,n}.Attributes.ID;
    StartTime = s.ExpoXData.Passes.Pass{1,n}.Attributes.StartTime;
    EndTime = s.ExpoXData.Passes.Pass{1,n}.Attributes.EndTime;
    BlockID = s.ExpoXData.Passes.Pass{1,n}.Attributes.BlockID;
    bid = str2num(BlockID);
    if bid == 0 || bid == 1 || bid == 2,
        m = m +1;
        block(m) = bid;
        slot(m) = str2num(SlotID);
        pass(m) = str2num(PassID);
        strT(m) = 0.1*str2num(StartTime);
        endT(m) = 0.1*str2num(EndTime);
        eind = size(s.ExpoXData.Passes.Pass{1,n}.Event,2);
        rid(m)  = str2num(s.ExpoXData.Passes.Pass{1,n}.Event{1,eind}.Attributes.RID);
        r_dat = s.ExpoXData.Passes.Pass{1,n}.Event{1,eind}.Attributes.Data;
        r = str2num(r_dat);
        rspt(m) = r(8)*1000 + strT(m);
    end
end
attri = cat(1,block,slot,rid,rspt,strT,endT,pass);
attri = attri';
end
%% code ends here

