function f_hkl=calc_f_hkl(f_j,h,k,l,x,y,z)
    f_hkl = 0;
    num_atoms = size(x,1);
    for j = 1:num_atoms
        f_hkl = f_hkl + f_j(:,j) .* exp( - 2*pi*1i * (h*x(j) + k*y(j)+ l*z(j) ));
    end
    %%
end