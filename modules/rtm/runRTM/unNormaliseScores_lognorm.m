%TPM 2020 - this function takes input NORMALISED cross-correlation peak heights
%(Z), and infers the raw XCF PH based on the probability distribution
%paramters pd.

%Z is a matrix of pattern number - by - trialled template structure.
%pds is a cell array of pattern number - by - trialled template structure,
%containing a 1-by-2 vector of [mu, sigma], lognormal parameters, in each
%cell.

%This is the inverse of NormaliseScores_lognorm

function [PH]=unNormaliseScores_lognorm(Z,pd)

PH=zeros(size(Z));
    for i=1:size(Z,1)
        for p=1:size(Z,2)
            
            %grab the peak height and the prob dist parameters
            Z0=Z(i,p);
            params=pd{i,p};

            %un=center the distribution
            Z0=(Z0*params(2))+params(1);
            
            %take the exp()
            exp_z=exp(Z0);
            PH(i,p)=exp_z;

        end
    end
end