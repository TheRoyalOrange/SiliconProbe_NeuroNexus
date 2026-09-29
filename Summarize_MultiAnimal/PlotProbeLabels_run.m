% PlotProbeLabels_run.m
%
% Description: Plots the channel/shank activity labels saved in ProbeInfo (by
%   ChannelShankLabelling_SponData.m / ChannelShankLabelling_StimData.m) for a list
%   of mice, using plotProbeLabels.m: one figure per label and mouse, with the shank
%   scores on top and the channel-score heatmap in the ProbeMaps layout below. The
%   probe for each mouse is the one whose ProbeInfo.Areas matches its region (as in
%   the summary scripts). Figures are not saved.
%
% Inputs:
%   User-set config:
%     mice (cell array of strings, 1 x nMice) - animal IDs (<animal>-ProbeInfo.mat)
%     region (cell array of strings, 1 x nMice) - brain area per mouse; must match exactly
%       one entry of that mouse's ProbeInfo.Areas (selects the probe)
%     labelnames (string array) - labels to plot; empty = every label found for that probe
%   E:\Roy\Processed Silicon Probe Data\ProbeInfo\<animal>-ProbeInfo.mat, per mouse
%     (fields: see plotProbeLabels.m)
%
% Outputs:
%   figs_labels (struct) - one field per mouse (valid field name of the animal ID), each the
%     figure handles returned by plotProbeLabels
%
% Dependencies: plotProbeLabels.m (same folder); ChannelShankLabelling_SponData.m /
%   ChannelShankLabelling_StimData.m (write the labels into ProbeInfo).

mice = {"20260423-p12"};
region = {"V1"};
labelnames = strings(1,0); %e.g. ["SponActivity","LightOnly_mixedIntensity"]; empty = all

assert(numel(region) == numel(mice), ...
    'region has %d entries but mice has %d (need one region per mouse)', numel(region), numel(mice))
%%
figs_labels = struct();
for m = 1:numel(mice)
    pinfo = load(fullfile(['E:\Roy\Processed Silicon Probe Data\ProbeInfo\' char(mice{m}) '-ProbeInfo.mat']),'ProbeInfo');
    prb = find(strcmp(pinfo.ProbeInfo.Areas, region{m}));
    assert(numel(prb) == 1, '%s: region "%s" matches %d probes in ProbeInfo.Areas (need exactly 1)', ...
        mice{m}, region{m}, numel(prb))
    figs_labels.(matlab.lang.makeValidName(char(mice{m}))) = plotProbeLabels(pinfo.ProbeInfo, prb, labelnames);
end
clear pinfo
