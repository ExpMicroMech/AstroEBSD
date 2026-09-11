function pattern = normaliseForImwrite(pattern)
    pattern(pattern<-3.3)=-3.3;
    pattern(pattern>3.3) = 3.3;
    pattern = rescale(pattern,0,1);
end