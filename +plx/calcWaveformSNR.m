function snr = calcWaveformSNR(waveforms)
% calculate SNR from Kelly et al 2007

n = size(waveforms,1);
mw = nanmean(waveforms);
peak   = nanmax(mw);
trough = nanmin(mw);
amp    = peak - trough;
r  = waveforms - repmat(mw, n,1);
sigma  = nanstd(r(:));
snr = amp/(2*sigma);

