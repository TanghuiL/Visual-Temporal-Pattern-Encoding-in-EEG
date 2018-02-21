function gfp = computeGFP(wave,tmw,nch)
% ------------------------------------------------------------------------
% Compute Global Field Power (gfp) from epoched EEG waves
% Input (by calling the readERP function):
% wave         - waveform data in 2D (channels x timepoints) array
% tmw          - epoch window [ms] of ERP waveform (e.g. -100:2:500 ms)
% Output:
% gfp          - global field power
%
% copyright (c) Huizhen Tang, e-mail: eagtang2007@gmail.com, Nov-7-2017
% ------------------------------------------------------------------------
gfp = [];
% %% Specify electrodes for gfp computation
    if nch == 0;
        chls = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
            21 22 23 24 25 26 27 28 29 30 31]; % chls specifies the electrodes included in the calculation
    else disp('Invalid input of nch!')
    end
wave = wave(chls, :); % exclude the unselected electrode from wave
% computing gfp for each time point and iterate through all the time points
for t = 1:length(tmw)
    aveu = mean(wave(:,t)); % aveu is the mean voltage at time point tmw(t)
    difu = []; 
    for n = 1:length(chls)
        difu(n) = (wave(n,t)-aveu)^2; 
    end
    gfp(t) = sqrt(sum(difu)/length(chls));
end

%% save the return gfp values
% % select folder for saving the file
% folder_dir = uigetdir;
% prompt = {'Enter the file name for the return gfp values'};
% dlg_title = 'file name';
% num_lines = 1;
% defAns = {'PR-90d-stdO'};
% answer = inputdlg(prompt,dlg_title,num_lines,defAns);
% % Abort if the user clicks 'Cancel'.
% if isempty(answer), disp('Aborted.');
%     return;
% end
% fname = answer{1};
% if ~status  %%%Handle empty value returned for unsuccessful conversion
%     msgbox('Invalid Input','Error in Parameter settings','error');
% end
% save([folder_dir '\' fname], 'gfp','wave') % save gfp to the selected folder

% %% plot waves and gfp for visual inspection
% figure
% plot(tmw,gfp,'b','LineWidth',4)
% legend('gfp','Location','northeast')
% hold on
% plot(tmw,wave(chls,:),'-m','LineWidth',0.5)
% plot(tmw,gfp,'b','LineWidth',4)
% title([fname '-gfp'])
% xlabel('Time (ms)');
% ylabel('Amplitude (µV)');
% hold off
% fig = gcf;
% fig.PaperUnits = 'inchs';
% fig.PaperPosition = [0 0 8 4];
% saveas (gcf,([fname '.png']));
% saveas (gcf,([fname '.fig']));

end

%% Code ends here.