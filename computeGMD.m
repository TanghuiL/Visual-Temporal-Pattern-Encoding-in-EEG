function gmd = computeGMD(type,nch)
% ------------------------------------------------------------------------
% Compute Global Map DISSimilarity (GMD/DISS) from epoched EEG waves
% Input (by calling the computeGFP function):
% type         - comparison between two time points or between conditions
%                 type = 1, between two time points t1 and t2
%                 type = 2, between two conditions c1 and c2
% wave         - waveform data in 2D (channels x timepoints) array
% tmw          - epoch window [ms] of ERP waveform (e.g. -100:2:500 ms)
% gfp          - global field power of each time point
% nch          - selected electrodes for gfp calculation(ie.EOG excluded)

% Output:
% gmd          - global map dissimilarity index
%
% gmd can be computed to describe the difference between two different time
% points within one condition or the same time point between two different
% conditions
%
% copyright (c) Huizhen Tang, e-mail: eagtang2007@gmail.com, Nov-7-2017
% ------------------------------------------------------------------------
%% gfp computation depending on the type selected
if type == 1 %type 1 gives gmd between two conditions at the same time point.
    
    % call readERP function to convert text (.dat) waveform to 2D array.
    [wave1 tmw1] = readERP; %Call the "readERP" function to convert dt file to 2D array into wave1
    [wave2 tmw2] = readERP; %Call again the function to convert dt file to 2D array into wave2
    % calculate gfp for each ERP
    nch = 0;
    gfp1 = computeGFP(wave1,tmw1,nch);
    gfp2 = computeGFP(wave2,tmw2,nch);
    if nch == 0;
        chls = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
            21 22 23 24 25 26 27 28 29 30 31]; % chls specifies the electrodes included in the calculation
    else disp('Invalid input of nch!')
    end
    
    gmd = []; aveu1 = []; aveu2 = []; u = []; v = [];
    for t = 1:length(tmw1)
        aveu1(t) = mean(wave1(:,t));   
        aveu2(t) = mean(wave2(:,t));
    end
    for t = 1:length(tmw1)
        difuv = [];
        for n = 1:length(chls)
            difuv(n) = ((wave1(n,t) - aveu1(t))/gfp1(t) - (wave2(n,t) - aveu2(t))/gfp2(t))^2;
        end
        gmd(t) = sqrt(sum(difuv)/n);
    end
    
elseif type == 2 % gmd between two time points in the same condition
    
    % call readERP function to convert text (.dat) waveform to 2D array.
    [wave tmw] = readERP; %The readERP function returns two arrays, wave & tmw
    
    %  call function gfp
    nch = 0;
    gfp = computeGFP(wave,tmw,nch); %the default value for nch should be 0
    if nch == 0;
        chls = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
            21 22 23 24 25 26 27 28 29 30 31]; % chls specifies the electrodes included in the calculation
    else disp('Invalid input of nch!')
    end
    
    gmd = []; aveu = [];
    for t = 1:length(tmw)
        aveu(t) = mean(wave(:,t));
    end
    for t1 = 1:length(tmw)
        for t2 = 1:length(tmw)
            difuv = [];
            for n = 1:length(chls)
                difuv(n) = ((wave(n,t1) - aveu(t1))/gfp(t1) - (wave(n,t2) - aveu(t2))/gfp(t2))^2;
            end
            gmd(t1,t2) = sqrt(sum(difuv)/n);
        end
    end
else disp('Invalid type!')
    return
end
end

%% Code ends here.