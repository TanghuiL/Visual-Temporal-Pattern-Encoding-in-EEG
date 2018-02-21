% ------------------------------------------------------------------------
% Compute Global Map DISSimilarity (GMD/DISS) from epoched EEG waves
% For multiple pair-comparions for VisA study.
%
% gmd are computed to describe the difference between two different time
% points within one condition or the same time point between two different
% conditions
%
% copyright (c) Huizhen Tang, e-mail: eagtang2007@gmail.com, Feb-21-2017
% ------------------------------------------------------------------------

clear all 

%% calculate 
gmd = computeGMD(type,nch);
