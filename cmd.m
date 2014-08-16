

import pdsa.*

load data/spikes1.mat
spikes1 = spikes;
load data/spikes2.mat
spikes2 = spikes;
clear spikes

[spikeStructs, spikesNew] = pdsa.concatenateSpikes(spikes1, spikes2)
