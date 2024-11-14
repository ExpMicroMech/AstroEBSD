function [Library_G,Library_PH,pd_params,lib]=fLibraryTest(template_library,library_G,FFTData_Ref,SettingsXCF,XCF_data_fill,usepar)
%Test the library for template matches and return the best orientation

library_size=size(library_G,3);

RegOut1=zeros(library_size,4);
roisize=SettingsXCF.roisize;
mesh=SettingsXCF.mesh;

if usepar == 1 %parallel search the library
    
    parfor n=1:library_size
        %Cross correlate and find the peak
        [RegOut1(n,:)] = fReg( FFTData_Ref,template_library(:,:,n),roisize,mesh,XCF_data_fill); %RegOut = [Xshift, Yshift, fullXCFheight, normXCFheight]
    end
    
else %do not paralle search
    
    for n=1:library_size
        %Cross correlate and find the peak
        [RegOut1(n,:)] = fReg( FFTData_Ref,template_library(:,:,n),roisize,mesh,XCF_data_fill); %RegOut = [Xshift, Yshift, fullXCFheight, normXCFheight]
    end
    
end

%find the max cross correlation
[Library_PH,max_index]=max(RegOut1(:,4));

%get some stats on the distribution
%[h,p]=chi2gof(log(RegOut1(:,4)));
%pd{1}=fitdist(RegOut1(:,4),'logLogistic');
pd_params(1)=mean(log(RegOut1(:,4)));
pd_params(2)=std(log(RegOut1(:,4)));

%find the orientation matrix
Library_G = library_G(:,:, max_index);

lib=RegOut1(:,4);
end