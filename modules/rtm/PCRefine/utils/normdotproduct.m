function ntd=normdotprod(A,B)
    ntd=sum(dot(A,B),'all')/(sqrt(sum(dot(A,A),'all'))*sqrt(sum(dot(A,A),'all')));
end
