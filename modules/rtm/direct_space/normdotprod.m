function ndp=normdotprod(A,B)
    ndp=sum(dot(A,B),'all')/(sqrt(sum(dot(A,A),'all'))*sqrt(sum(dot(B,B),'all')));
end