function [ m_plot,I] =  cleanFHKLDuplicates(m_plot,I,InputUser)
   

    if nargin ==3
        try
            numHkl = InputUser.numHKLPLot;
            numHklClean = ceil(numHkl/2);
        catch
            warning('more than 2 arguments specified, but none contains maximum plot bands')
        end
    else
        numHkl = 15;
        numHklClean = 8;
    end
    try
    m_plot = m_plot(1:numHkl);

    I = I(1:numHkl);
    % catch

    end
    % parallel = angle(m_plot,m_plot')
    indices = (1:numHkl);
    %%
    dropList = [];

    for m_index =1: numel(m_plot)
        maxIIndex = find(angle( m_plot(m_index), m_plot)==0);
        % maxIIndex(maxIIndex==m_index)=[];
        I_drop = find(I(maxIIndex) ~=max(I(maxIIndex)));
        % I(I_drop) =[];
        % m_plot()
        dropList = [dropList maxIIndex(I_drop)'];
        
            % disp('duplicate')

        % end
    end

    %%
    I(dropList) = [];
    m_plot(dropList) =[];
    %%
    if m_plot.CS.numSym>=24 %cubic :D
        try
        I = I(1:5);
        m_plot = m_plot(1:5);
        catch
            warning('m_plot allready small enough')
        end
    else
        try
            I = I(1:numHklClean);
            m_plot = m_plot(1:numHklClean);
        catch
            try
            I = I(1:10);
            m_plot = m_plot(1:10);
            catch
            warning('m_plot allready small enough')
            end
        end
    end
    

end