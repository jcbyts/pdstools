function dp=dprime(x,y)
    dp=(mean(x(y,:))-mean(x(~y,:)))./sqrt((var(x(y,:))+var(x(~y,:)))/2);