function [ m_plot,I] cleanFHKLDuplicates()


    numHkl =15;
    m_plot = m_plot(1:numHkl);
    I = I(1:numHkl);
    % parallel = angle(m_plot,m_plot')
    indices = (1:numHkl);
    %%
    dropList = [];
    for m_index = size(m_plot)
        maxIIndex = find(angle( m_plot(m_index), m_plot(indices~=m_index))==0);
        I_drop = find(I(maxIIndex) ~=max(I(maxIIndex)));
        % I(I_drop) =[];
        % m_plot()
        dropList = [dropList maxIIndex(I_drop)];
        
            % disp('duplicate')

        % end
    end
    I(dropList) = [];
    m_plot(dropList) =[];
    %%
    if m_plot.CS.numSym>=24 %cubic :D
        I = I(1:5);
        m_plot = m_plot(1:5);
    else
        I = I(1:10);
        m_plot = m_plot(1:10)
    end
    
    structConfig.IPlot = I;
end