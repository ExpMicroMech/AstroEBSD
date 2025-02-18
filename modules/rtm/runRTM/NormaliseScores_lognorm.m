%TPM 2020 - this function takes input cross-correlation peak heights
%(PH_Out), and ingers the Z-score based on the probability distribution
%paramters pd.

%PH_Out is a matrix of pattern number - by - trialled template structure.
%pds is a cell array of pattern number - by - trialled template structure,
%containing a 1-by-2 vector of [mu, sigma], lognormal parameters, in each
%cell.

function [Z]=NormaliseScores_lognorm(PH_Out,pd)

Z=zeros(size(PH_Out));
    for i=1:size(PH_Out,1)
        for p=1:size(PH_Out,2)
            
            %grab the peak height and the prob dist parameters
            PH=PH_Out(i,p);
            params=pd{i,p};

            %take the exp() and then mean centre / normalise
            log_ph=log(PH);
            log_ph_norm=(log_ph-params(1))./params(2);
            Z(i,p)=log_ph_norm;

        end
    end
end