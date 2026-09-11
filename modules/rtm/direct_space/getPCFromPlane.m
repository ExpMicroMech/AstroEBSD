function PC = getPCFromPlane(x,y,PCX_mod,PCY_mod,PCZ_mod)
% PCX_mod,PCY_mod,PCZ_mod
PC = zeros(numel(x),3);
PC(1) = x*PCX_mod(1)+y*PCX_mod(2)+PCX_mod(3);
PC(2) = x*PCY_mod(1)+y*PCY_mod(2)+PCY_mod(3);
PC(3) = x*PCZ_mod(1)+y*PCZ_mod(2)+PCZ_mod(3);
end