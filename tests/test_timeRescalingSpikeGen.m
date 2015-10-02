%% test_timeRescalingSpikeGen
% use time rescaling theorem to transform inter spike intervals from a
% gamma process into an over or under dispersed inhomogeneous poisson
% process

% set up parameters
bs       = 1e-3; % bin size (seconds)
time     = linspace(0,10, 10/bs);
lambda   = (10*sin(2*pi*time)+15); % spike rate
nRepeats = 100;
% shape parameter of gamma distribution
% 1 > super poisson
% 1 = poisson
% 1 < sub poisson
shapes=[.02 1 5]; 
labels={'superpoisson', 'poisson', 'subpoisson'};

for kShape=1:3
    
    spks = zeros(nRepeats, numel(time));
    for k = 1:nRepeats
        spks(k,:) = genSpikesRenewalProcess(lambda, shapes(kShape), bs);
    end

    [i,j] = find(spks);
    figure(kShape); clf
    subplot(3,1,1)
    plot(time(j), i, '.')
    xlabel('time')
    ylabel('trial')
    title(labels{kShape})
    
    subplot(3,1,2)
    smwin = 40;
    plot(time,smooth(mean(spks),smwin)/bs, 'k');hold on
    plot(time,smooth(var(spks),smwin)/bs, '-');hold on
    plot(time,lambda, 'r', 'Linewidth', 2)
    ylim([0 100])
    legend({'psth', 'pstv', 'true rate'})
    xlabel('time')
    ylabel('spike rate')
    title(sprintf('%d ms smoothing window', smwin*bs*1e3))
    
   
    subplot(3,1,3)
    plot(time,smooth(var(spks),smwin)./smooth(mean(spks),smwin))
    hold on
    plot(xlim, [1 1], 'k--')
    ylabel('fano factor')
    title('Variance Mean Ratio')
    xlabel('time')
end