% plotProbeLabels.m
%
% Description: Redraws the channel/shank activity-label plots of
%   ChannelShankLabelling_SponData.m and ChannelShankLabelling_StimData.m from
%   what those scripts save in ProbeInfo, so they can be made without rerunning
%   the labelling. One figure per label: top = shank score per shank (line with
%   markers, one point per ProbeMaps column), bottom = heatmap of the channel
%   scores in the ProbeMaps layout (columns = shanks, rows = depth, row 1 at the
%   top), with each cell labelled by its channel ID and Ch_Remove channels shown
%   black with an 'X'. Same style as the labelling scripts: orange colour map for
%   spontaneous labels (name contains "Spon"), green for stimulus-driven ones.
%   Returns the figure handles and saves nothing (meant to be reused later next to
%   plots of the summary csv / tableres).
%
% Inputs:
%   ProbeInfo (struct) - an animal's ProbeInfo (<animal>-ProbeInfo.mat). Fields used:
%     .ProbeMaps (cell, 1 x nProbes) - rows x shanks raw channel IDs (0 = no site)
%     .ShankLabels (cell, 1 x nProbes, optional) - per probe a struct of shank labels,
%       each 1 x shanks (ProbeMaps column order)
%     .ChanLabels (cell, 1 x nProbes, optional) - per probe a struct of channel labels,
%       each rows x shanks (ProbeMaps layout; NaN = no score, e.g. removed channel)
%     .chLabels (cell, 1 x nProbes, optional) - older name for ChanLabels (used by earlier
%       versions of ChannelShankLabelling_SponData.m); read when a label isn't in ChanLabels
%     .Ch_Remove (cell, 1 x nProbes, optional) - raw channel IDs removed per probe
%     .Animal (char, optional), .Areas (cell, optional) - used in the titles
%   prb (double, scalar) - probe index (into ProbeMaps / the label cells)
%   labelnames (string array, optional) - labels to plot; default = every label found for
%     this probe in ShankLabels, ChanLabels and chLabels (in that order)
%   parent (TiledChartLayout, optional) - draw into tiles of this layout instead of making new
%     figures (used by plotLWIntensityCompare.m); each label gets a nested 4x1 layout there
%   tiles (double, 1 x numel(labelnames), required with parent) - tile number of each label in parent
%
% Outputs:
%   figs (graphics handle array) - one figure per label plotted (or, with parent, one nested
%     layout per label); empty if the probe has no labels
%
% Dependencies: ChannelShankLabelling_SponData.m / ChannelShankLabelling_StimData.m
%   (write ProbeInfo.ChanLabels / ShankLabels); OpenEphys_BaseAnalysis*.m (ProbeMaps);
%   OpenEphys_editProbeInfo_ChRemove.m (Ch_Remove).

function figs = plotProbeLabels(ProbeInfo, prb, labelnames, parent, tiles)

figs = gobjects(1,0);
useparent = nargin >= 4 && ~isempty(parent);
pmap = ProbeInfo.ProbeMaps{prb};  %rows x shanks raw channel IDs
[nrow,nshnk] = size(pmap);
removed = false(size(pmap));
if isfield(ProbeInfo,'Ch_Remove') && numel(ProbeInfo.Ch_Remove) >= prb
    removed = ismember(pmap, ProbeInfo.Ch_Remove{prb});
end

SL = getLabels(ProbeInfo,'ShankLabels',prb);
CL = getLabels(ProbeInfo,'ChanLabels',prb);
CLold = getLabels(ProbeInfo,'chLabels',prb);
found = unique([string(fieldnames(SL))' string(fieldnames(CL))' string(fieldnames(CLold))'],'stable');
if nargin < 3 || isempty(labelnames)
    labelnames = found;
else
    labelnames = string(labelnames);
end
if isempty(labelnames)
    fprintf('plotProbeLabels: probe %d has no channel/shank labels - nothing to plot\n', prb);
    return
end
if useparent
    assert(nargin >= 5 && numel(tiles) == numel(labelnames), 'plotProbeLabels: give one tile number per label with parent')
end

animal = '';
if isfield(ProbeInfo,'Animal'), animal = char(ProbeInfo.Animal); end
area = '';
if isfield(ProbeInfo,'Areas') && numel(ProbeInfo.Areas) >= prb, area = char(ProbeInfo.Areas{prb}); end

%colour maps of the labelling scripts (never reach black, so removed channels stay distinct)
orangemap = interp1([0 0.5 1],[1.00 0.96 0.90; 1.00 0.62 0.15; 0.92 0.38 0.00],linspace(0,1,256)); %spontaneous
greenmap = interp1([0 0.5 1],[0.95 0.98 0.93; 0.55 0.82 0.40; 0.18 0.60 0.22],linspace(0,1,256));  %stimulus-driven

for L = 1:numel(labelnames)
    lab = char(labelnames(L));

    %shank scores (1 x shanks)
    shk = [];
    if isfield(SL,lab)
        shk = double(SL.(lab)(:)');
        if numel(shk) ~= nshnk
            warning('plotProbeLabels: ShankLabels{%d}.%s has %d values but the probe has %d shanks - not plotted', prb, lab, numel(shk), nshnk)
            shk = [];
        end
    end
    %channel scores (ProbeMaps layout)
    chm = [];
    if isfield(CL,lab)
        chm = CL.(lab);
    elseif isfield(CLold,lab)
        chm = CLold.(lab);
    end
    if ~isempty(chm)
        if ~isequal(size(chm),size(pmap)) || ~isnumeric(chm)
            warning('plotProbeLabels: channel label %s of probe %d is not a numeric %dx%d map - not plotted', lab, prb, nrow, nshnk)
            chm = [];
        else
            chm = double(chm);
            chm(removed) = NaN;
        end
    end
    if isempty(shk) && isempty(chm)
        warning('plotProbeLabels: label "%s" not found for probe %d', lab, prb)
        continue
    end

    if contains(lab,'Spon','IgnoreCase',true)
        cmap = orangemap;
    else
        cmap = greenmap;
    end
    if useparent
        %nested 4x1 layout in the given tile of the parent layout
        tl = tiledlayout(parent,4,1,'TileSpacing','compact');
        tl.Layout.Tile = tiles(L);
        title(tl,lab,'Interpreter','none');
        figs(end+1) = tl; %#ok<AGROW>
    else
        ttl = sprintf('%s - probe %d (%s): %s', animal, prb, area, lab);
        figs(end+1) = figure('Name',ttl,'Color','w'); %#ok<AGROW>
        tl = tiledlayout(figs(end),4,1,'TileSpacing','compact');
        title(tl,ttl,'Interpreter','none');
    end

    %top: shank scores, one point per ProbeMaps column (sits above its heatmap column)
    ax1 = nexttile(tl,1);
    if ~isempty(shk)
        plot(ax1,1:nshnk,shk,'-o','LineWidth',2,'MarkerSize',8,'Color',cmap(end,:),'MarkerFaceColor',cmap(end,:));
    else
        text(ax1,0.5,0.5,'no shank scores for this label','Units','normalized','HorizontalAlignment','center');
    end
    xlim(ax1,[0.5 nshnk+0.5]);
    xticks(ax1,1:nshnk); xticklabels(ax1,{});
    ylabel(ax1,{'Shank score','(mean of channels)'});
    grid(ax1,'on'); ax1.GridAlpha = 0.15; box(ax1,'off');

    %bottom: channel heatmap (NaN cells are transparent over a black axes background)
    ax2 = nexttile(tl,2,[3 1]);
    if isempty(chm)
        chm = NaN(nrow,nshnk); %layout only
        title(ax2,'no channel scores for this label','FontWeight','normal');
    end
    imagesc(ax2,1:nshnk,1:nrow,chm,'AlphaData',~isnan(chm));
    ax2.Color = 'k';
    colormap(ax2,cmap);
    if any(~isnan(chm),'all')
        cb = colorbar(ax2);
        cb.Layout.Tile = 'east'; %outside both tiles, so the shank columns stay aligned
        cb.Label.String = [lab ' (channel score, a.u.)'];
        cb.Label.Interpreter = 'none';
    end
    xticks(ax2,1:nshnk); yticks(ax2,1:nrow);
    xlabel(ax2,'Shank (ProbeMaps column)'); ylabel(ax2,'ProbeMaps row');
    linkaxes([ax1 ax2],'x');

    %cell labels: channel ID, or 'X' for removed channels; white text on black cells, black on coloured
    for r = 1:nrow
        for c = 1:nshnk
            if pmap(r,c) == 0
                continue %no site
            elseif removed(r,c)
                txt = 'X';
            else
                txt = num2str(pmap(r,c));
            end
            if isnan(chm(r,c))
                txtcol = 'w';
            else
                txtcol = 'k';
            end
            text(ax2,c,r,txt,'HorizontalAlignment','center','Color',txtcol,'FontSize',8);
        end
    end
end
end

function S = getLabels(ProbeInfo, fieldname, prb)
%struct of labels for probe prb from ProbeInfo.(fieldname), or an empty struct if there is none
S = struct();
if isfield(ProbeInfo,fieldname) && iscell(ProbeInfo.(fieldname)) && numel(ProbeInfo.(fieldname)) >= prb ...
        && isstruct(ProbeInfo.(fieldname){prb})
    S = ProbeInfo.(fieldname){prb};
end
end
