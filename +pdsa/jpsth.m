function S = jpsth(s1,s2,lag)
% calculate jpsth (Brody, 1999)
% jpsth(s1,s2)
% Inputs
%	s1 [nTrials x nBins] binned spike counts for neuron 1
%	s2 [nTrials x nBins] binned spike counts for neuron 2
	if ~exist('lag', 'var')
		lag=50;
	end
	m1=mean(s1);
	m2=mean(s2);
	sd1=std(s1);
	sd2=std(s2);

	jpnorm= (rawJPSTH(s1,s2)-(m1(:)*m2(:)')) ./ (sd1(:)*sd2(:)');
	[covgrm,shigh,slow, ccraw]=covariogram(s1,s2,lag);
	speak=sigspan(covgrm,shigh);
	strough=sigspan(-covgrm,-slow);

	S=struct('psth1',m1,'psth2',m2,...
		'njpsth', jpnorm, ...
		'xchist', xcorrhist(jpnorm,lag), ...
        'ccraw', ccraw, ...
		'covgrm', covgrm,...
		'slow',slow, 'shigh', shigh, ...
		'speak', speak, 'strough', strough);


function rawjp=rawJPSTH(s1,s2)
	s=size(s1,2);
	assert(s==size(s2,2), 'Binned window for each neuron must be the same size')
	rawjp=zeros(s);
	for ii=1:s
		for jj=1:s
			rawjp(ii,jj) = mean(s1(:,ii).*s2(:,jj));
		end
	end

function xch=xcorrhist(jp,lag)
	if nargin < 2, lag=50; end
	xch=zeros(1,lag-lag+1);
	for k=-lag:lag
		xch(k+lag+1) = (sum(diag(jp,k)))/(length(jp)-abs(k));
	end

function [covgrm, shigh, slow, cc] = covariogram(s1,s2,lag)
	if nargin<3, lag=50; end
	nTrials=size(s1,1);
	nBins=size(s1,2);
	n=2*lag+1;
	v1v2=zeros(1,n);
	m1v2=zeros(1,n);
	v1m2=zeros(1,n);		
	cc=zeros(1,n);
	ccshuff=zeros(1,n);

	m1=mean(s1);
	m2=mean(s2);
	sd1=std(s1);
	sd2=std(s2);

	for k = 1:n
		currLag=k-lag-1;
		if currLag < 0
			jvec=1-currLag:nBins;
		else
			jvec=1:nBins-currLag;
		end

		for j = jvec
			cc(k)=cc(k)+mean(s1(:,j).*s2(:,j+currLag));
			ccshuff(k)=ccshuff(k) + m1(j) * m2(j+currLag);
			v1v2(k)=v1v2(k) + sd1(j).^2 * sd2(j+currLag).^2;
			m1v2(k)=m1v2(k) + m1(j).^2 * sd2(j+currLag).^2;
			v1m2(k)=v1m2(k) + sd1(j).^2 * m2(j+currLag).^2;
		end
	end

	covgrm=cc-ccshuff;
	sigma=sqrt((v1v2 + m1v2 + v1m2)/nTrials);
	shigh=2*sigma;
	slow=-2*sigma;

function [spanE] = sigspan(x, sig)
	% TODO: reimplement with spanlocs.m
	if numel(sig)==1
		sig = repmat(sig,1,numel(x));
	end
		
	d = x > sig;

    dd = [d(2:end) 0] - d(:)';
        
    bix = find(dd==1) + 1; % beginning indices
    eix = find(dd==-1); % ending indices

    if ~isempty(eix) && length(bix) < length(eix) 		 
        bix = [1 bix]; 		 
    end

    % identify the longest span(s)
	maxIndices = (eix-bix)==max(eix-bix);
	
    %find the indices of the longest spans
    spanE = [];
    theIndices = find(maxIndices==1);
    for i=1:length(theIndices)
        theIndex = theIndices(i);
        spanE = horzcat(spanE,[bix(theIndex) eix(theIndex)]); %#ok<AGROW>
    end