function extend = findMaxExtend(nVec)
    excludes = nVec.uvw < 1.e-2  ;
    % [u v w] < 1e-2
    uvw = [nVec.uvw];
    uvw(excludes) =1;
    u = uvw(1);
    v = uvw(2);
    w = uvw(3);
    results = [ abs(linsolve(u,1)),abs(linsolve(v,1)),abs(linsolve(w,1)) ];
    results(excludes) = 0;
    % c = abs(linsolve(w,1));

    extend =(max(results(:)));
   %%
end