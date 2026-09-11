function ncc =  nxCorrelation(A,B)
%     A = rescale(A)
%     B = rescale(B);
% % %     meanA = mean(A(:));
% % %     meanB = mean(B(:)); % test if we get some wird numeric error by subtracting with a small numer
% % %     
% % %     if meanA<1e-5
% % %         meanA = 0;
% % %     end
% % %     if meanB <1e-5
% % %         meanB = 0;
% % %     end
% corrcoef(A,B)
%%
    normA = sqrt(sum( (A(:) - mean(A(:)) ).^2 ));
    normB = sqrt(sum( (B(:) - mean(B(:)) ).^2 ));
    ncc = sum( (A(:) - mean(A(:))) .*( B(:) -mean(B(:)) )) ...
                /( normA * normB);
    % disp(class(A))
    % imwrite(rescale(A),'A.png')
    % imwrite(rescale(B),'B.png')
% disp(ncc);
    %%
end