function binnedData = binData(dataVector, binSize)
% binnedData = binData(dataVector, binSize)

if binSize == 1
    binnedData = dataVector;
    return
end

a = size(dataVector);
if a(2)>a(1)
    dataVector = dataVector';
    a = fliplr(a);
end

n=ceil(a(1)/binSize);
binnedData=zeros(n,1);
for t=1:(n-1)
    binnedData(t)=mean(dataVector((t-1)*binSize+(1:binSize)));
end