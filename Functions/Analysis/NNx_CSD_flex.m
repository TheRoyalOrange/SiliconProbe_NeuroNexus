% NNx_CSD_flex.m
%
% Description: Inverse (delta) current source density along each shank of
%   one probe, for every trial and time point at once. Same equation as
%   CSDallshank.m (Pettersen et al. 2006, J Neurosci Methods 154:116-33),
%   but contact depths come from ProbeInfo.ProbeMaps, so shanks with
%   Ch_Remove channels missing are handled at their real depths instead of
%   being closed up. Each contact's slab thickness h (the stretch of shank
%   one contact stands for) is half the gap to the contact above plus half
%   the gap to the contact below; the open side of an end contact counts
%   as half the nominal spacing. On a complete shank this gives h = spacing
%   for every contact, i.e. exactly CSDallshank.m. No plotting, no saving.
%
%   Usage:
%     [csd, csdinfo] = NNx_CSD_flex(data, ProbeInfo, prb)                    %all probe channels
%     [csd, csdinfo] = NNx_CSD_flex(data, ProbeInfo, prb, 'ChRemove', true)  %Ch_Remove already dropped
%
% Inputs:
%   data (double/single, trials x time x channels) - LFP in microvolts.
%     Channel order INFERRED to match the pipeline:
%       'ChRemove' false: data(:,:,k) = ProbeInfo.ProbeIds{prb}(k)
%       'ChRemove' true:  same list with ProbeInfo.Ch_Remove{prb} dropped
%                         (= chansuse in NNx_ITPC_PhaseReset_singleanimal.m)
%     A single trial must still be 1 x time x channels.
%   ProbeInfo (struct) - from NNx_BaseAnalysis.m; uses
%     .ProbeIds{prb} (raw channel IDs on the probe), .Ch_Remove{prb} (raw
%     channel IDs to exclude; may be missing for this probe), .ProbeMaps{prb}
%     (2D [depth rows x shanks] matrix of raw channel IDs). INFERRED: rows
%     are equally spaced contacts in physical order along the shank.
%   prb (int) - probe index into the ProbeInfo cell fields
%   Name-value options:
%     'ChRemove' (logical, default false) - false: data holds every probe
%       channel and every shank is treated as complete. true: Ch_Remove
%       channels are not in data; depths of the rest come from ProbeMaps.
%     'spacing' (double, default 100) - distance between adjacent ProbeMaps
%       rows, um
%     'radius' (double, default 177) - radius of the current disc around
%       each contact, um (default = the CSDallshank call in
%       NNx_BaseAnalysis.m)
%     'conductivity' (double, default 0.3) - extracellular conductivity, S/m
%     'minchans' (int, default 3) - shanks with fewer channels in data get
%       NaN CSD
%
% Outputs:
%   csd (double, same size and channel order as data) - CSD in uA/mm^3.
%     Sign as in CSDallshank.m (positive = source, negative = sink).
%   csdinfo (struct) - fields:
%     .chan_ids (nch x 1) - raw channel ID of each data channel
%     .shank, .depths_um, .h_um (nch x 1) - shank (ProbeMaps column), depth
%       (row x spacing, um) and slab thickness (um) of each data channel
%     .nchans_shank (1 x nshank) - channels in data per shank
%     .nan_shanks (int vector) - shanks left NaN (fewer than minchans)
%     .ChRemove, .spacing_um, .radius_um, .conductivity, .minchans - options used
%     .units (char) - 'uA/mm^3'
%
% Dependencies: none (ProbeInfo from NNx_BaseAnalysis.m,
%   Ch_Remove from NNx_editProbeInfo_ChRemove.m).

function [csd, csdinfo] = NNx_CSD_flex(data, ProbeInfo, prb, varargin)

p = inputParser;
addRequired(p, 'data', @(x) isnumeric(x) && ndims(x) <= 3);
addRequired(p, 'ProbeInfo', @isstruct);
addRequired(p, 'prb', @(x) isscalar(x) && x >= 1);
addParameter(p, 'ChRemove', false, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
addParameter(p, 'spacing', 100, @(x) isscalar(x) && x > 0);
addParameter(p, 'radius', 177, @(x) isscalar(x) && x > 0);
addParameter(p, 'conductivity', 0.3, @(x) isscalar(x) && x > 0);
addParameter(p, 'minchans', 3, @(x) isscalar(x) && x >= 1);
parse(p, data, ProbeInfo, prb, varargin{:});
opt = p.Results;
chremove = logical(opt.ChRemove);

%raw channel IDs in data order
probeids = ProbeInfo.ProbeIds{prb}(:);
removed = [];
if numel(ProbeInfo.Ch_Remove) >= prb
    removed = ProbeInfo.Ch_Remove{prb}(:);
end
if chremove
    chan_ids = probeids(~ismember(probeids, removed));
else
    chan_ids = probeids;
end

%data must have the channel count the ChRemove option implies
nch = size(data,3);
if nch ~= numel(chan_ids)
    if ~chremove
        error('NNx_CSD_flex:channelCount', ['data has %d channels, but probe %d has %d channels in ProbeInfo.ProbeIds ' ...
            '(''ChRemove'', false expects all of them). If the %d Ch_Remove channel(s) were already dropped ' ...
            'from data, call NNx_CSD_flex with ''ChRemove'', true.'], nch, prb, numel(probeids), numel(removed));
    elseif nch == numel(probeids)
        error('NNx_CSD_flex:channelCount', ['data has %d channels, the full probe %d channel count, but ''ChRemove'', true ' ...
            'expects %d (Ch_Remove dropped). If data still holds every channel, call NNx_CSD_flex with ''ChRemove'', false.'], ...
            nch, prb, numel(chan_ids));
    else
        error('NNx_CSD_flex:channelCount', ['data has %d channels, but probe %d has %d channels after dropping ' ...
            'its %d Ch_Remove channel(s) (''ChRemove'', true).'], nch, prb, numel(chan_ids), numel(removed));
    end
end

probemap = ProbeInfo.ProbeMaps{prb}; %[depth rows x shanks]
if ~all(ismember(chan_ids, probemap(:)))
    error('NNx_CSD_flex:probeMap', 'Channel(s) %s of probe %d are not in ProbeInfo.ProbeMaps{%d}.', ...
        mat2str(chan_ids(~ismember(chan_ids, probemap(:)))'), prb, prb);
end
nshank = size(probemap,2);
ntr = size(data,1);
nt = size(data,2);

r = opt.radius*1e-6; %m
csd = nan(size(data));
shank = nan(nch,1);
depths_um = nan(nch,1);
h_um = nan(nch,1);
nchans_shank = zeros(1,nshank);
nan_shanks = [];
for q = 1:nshank
    rows = find(ismember(probemap(:,q), chan_ids)); %rows of this shank present in data, in depth order
    [~, loc] = ismember(probemap(rows,q), chan_ids); %their positions in data
    nq = numel(rows);
    nchans_shank(q) = nq;
    shank(loc) = q;
    z = rows(:)*opt.spacing; %depth, um
    depths_um(loc) = z;
    gaps = diff(z);
    hq = ([opt.spacing; gaps] + [gaps; opt.spacing])/2; %slab: half gap above + half gap below, open side = spacing/2
    h_um(loc) = hq;
    if nq < opt.minchans
        nan_shanks(end+1) = q; %#ok<AGROW>
        continue
    end

    %V = F*C, F(j,i) = h_i/(2*sigma) * (sqrt((z_j-z_i)^2 + r^2) - |z_j-z_i|)
    zm = z*1e-6; %m
    dz = zm - zm'; %(j,i) = z_j - z_i
    F = (sqrt(dz.^2 + r^2) - abs(dz)) .* (hq'*1e-6/(2*opt.conductivity));
    V = reshape(double(data(:,:,loc)), ntr*nt, nq)' * 1e-6; %channels x (trials*time), uV -> V
    C = (F \ V) * 1e-3; %A/m^3 -> uA/mm^3
    csd(:,:,loc) = reshape(C', ntr, nt, nq);
end

csdinfo = struct();
csdinfo.chan_ids = chan_ids;
csdinfo.shank = shank;
csdinfo.depths_um = depths_um;
csdinfo.h_um = h_um;
csdinfo.nchans_shank = nchans_shank;
csdinfo.nan_shanks = nan_shanks;
csdinfo.ChRemove = chremove;
csdinfo.spacing_um = opt.spacing;
csdinfo.radius_um = opt.radius;
csdinfo.conductivity = opt.conductivity;
csdinfo.minchans = opt.minchans;
csdinfo.units = 'uA/mm^3';

end
