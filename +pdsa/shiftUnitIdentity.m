function spikeStructs = shiftUnitIdentity(spikeStructs, structId, spikeId, newSpikeId)
% Shift spike Unit id in a spikes struct
% spikeStructs = shiftUnitIdentity(spikeStructs, structId, spikeId, newSpikeId)

idx1 = spikeStructs{structId}.id == spikeId;
if isnan(newSpikeId)
    spikeStructs{structId}.time(idx1) = [];
    spikeStructs{structId}.id(idx1) = [];
    spikeStructs{structId}.waveform(idx1,:) = [];
    spikeStructs{structId}.snr(spikeId) = nan;
    spikeStructs{structId}.channel(spikeId) = nan;
    return
end

idx2 = spikeStructs{structId}.id == newSpikeId;
snr1 = spikeStructs{structId}.snr(spikeId);

ch1 = spikeStructs{structId}.channel(spikeId);
if newSpikeId > numel(spikeStructs{structId}.snr)
    snr2 = nan;
    ch2 = nan;
else
    snr2 = spikeStructs{structId}.snr(newSpikeId);
    ch2 = spikeStructs{structId}.channel(newSpikeId);
end

% change them over
spikeStructs{structId}.id(idx1) = newSpikeId;
spikeStructs{structId}.id(idx2) = spikeId;
spikeStructs{structId}.snr(spikeId)     = snr2;
spikeStructs{structId}.snr(newSpikeId)  = snr1;
spikeStructs{structId}.channel(spikeId) = ch2;
spikeStructs{structId}.channel(newSpikeId) = ch1;
