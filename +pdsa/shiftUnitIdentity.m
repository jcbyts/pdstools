function spikeStructs = shiftUnitIdentity(spikeStructs, structId, spikeId, newSpikeId);

idx1 = spikeStructs{structId}.id == spikeId;
idx2 = spikeStructs{structId}.id == newSpikeId;
snr1 = spikeStructs{structId}.snr(spikeId);
snr2 = spikeStructs{structId}.snr(newSpikeId);
ch1 = spikeStructs{structId}.channel(spikeId);
ch2 = spikeStructs{structId}.channel(newSpikeId);

% change them over
spikeStructs{structId}.id(idx1) = newSpikeId;
spikeStructs{structId}.id(idx2) = spikeId;
spikeStructs{structId}.snr(spikeId)     = snr2;
spikeStructs{structId}.snr(newSpikeId)  = snr1;
spikeStructs{structId}.channel(spikeId) = ch2;
spikeStructs{structId}.channel(newSpikeId) = ch1;