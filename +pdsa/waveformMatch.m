function [a, d1, d2] = waveformMatch(x,y)
% Waveform matching: from Tolias et al 2007
% Use a multistep procedure to compare waveforms for each neuron.
% I replaced their amplitude measure with the pearson's correlation 
% coefficient. 
% For average waveforms x, y
% 1) scale x to minimize the least squares distance between x and y
% 
% 		a(x,y) = argmin a  ||ax - y ||^2                     (7)
% 2) compute euclidean distance between scaled waveforms
% 		d1(x,y) = sum(  ||a(x, y)x - y|| / ||y|| 			 (8)
% 3) 

nx = size(x,2);
ny = size(y,2);

a  = zeros(nx, ny);
d1 = zeros(nx, ny);
d2 = zeros(nx, ny);
for ii = 1:nx
	x(:,ii) = x(:,ii) - mean(x(:,ii)); 
	for jj = 1:ny
		y(:,jj) = y(:,jj) - mean(y(:,jj)); 
		a(ii,jj) = (x(:,ii)'*x(:,ii))\(x(:,ii)'*y(:,jj));
		d1(ii,jj) = norm(a(ii,jj)*x(:,ii) - y(:,jj), 2)/norm(y(:,jj), 2);
		d2(ii,jj) = (x(:,ii)'*y(:,jj)) ./ ( sqrt(sum(x(:,ii).^2)) * sqrt(sum(y(:,jj).^2)));
	end
end

% d2 = (x'*y) ./ bsxfun(@times, sqrt(sum(x.^2))', sqrt(sum(y.^2)))
