function ndp=dotproduct_optimization(pcx,pcy,pcz,phi1,PHI,phi2,screen_int,Settings_Cor,Pattern_in)
    compare_pat=generate_pattern(screen_int,Settings_Cor,pcx,pcy,pcz,phi1,PHI,phi2)
    compare_pat=EBSP_BGCor(compare_pat)
    ndp=normdotprod(normalize(PC_in),normalize(compare_pat))
end