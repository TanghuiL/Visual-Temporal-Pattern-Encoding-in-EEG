function fig = PlotPSTH(m, sem, t, timeUnit, spikeUnit, fig, plotArgs)
% function fig = PlotPSTH(m, sem, t, timeUnit, spikeUnit, fig, plotArgs)
%
% Plot a post-stimulus-time histogram with the specified mean and standard
% error (sem).  t specifies the time at which each bin begins.
% [optional:]
% timeUnit: is printed below the x-axis as the unit of time.  Should name the units of t [default: ''].
% spikeUnit: is printed as the y-axis unit.  Should name the units of the m, sem [default: ''].
% fig: handle of existing figure to plot to. [default: make a new one].
%   Fig can also be an array [fig ax], in which case the two elements are a
%   figure handle and an axis handle, and the plot is on those axes.
% plotArgs are passed to plot [default: 'k'].
% Returns fig, or the newly-created fig.
%
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetSpikeTimes, GetAnalog,
% GetSpikeCounts, GetPSTH, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor.
%
%   Author:      Jim Muller
%   Last updated:  2005-01-11
%   E-mail:      jim@monkeybiz.stanford.edu


if ~exist('timeUnit'), timeUnit = ''; end
if ~exist('spikeUnit'), spikeUnit = ''; end
if ~exist('fig')
    fig=figure; 
else 
    if ~isscalar(fig)
        figure(fig(1));
        axes(fig(2));
    else
        figure(fig);
    end
end
if ~exist('plotArgs'), plotArgs = 'k'; end

hold on
set(gca,'FontSize', 16, 'FontWeight', 'Bold')
plot(t, m, plotArgs, 'LineWidth', 2);
[errx, erry] = MakeErrorBars(t, m, sem);
plot(errx, erry, plotArgs, 'LineWidth', 2);
xlabel(['Time (' timeUnit ')'])
ylabel(['Response (' spikeUnit ')'])
