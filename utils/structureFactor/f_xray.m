function f=f_xray(a_i,b_i,c, theta,lambda)
%%
    k = 2* pi / lambda;
    Q = 2* k * sin(theta);
    Q = reshape(Q, 1,[]); % align vector, such thta matrix mult works
    f=0;
    for i =1:4
        f= f + a_i(:,i) .* exp( -b_i(:,i) .* (Q/ (4*pi)).^2 ) + c(:);
    % f
    end
    %%
end