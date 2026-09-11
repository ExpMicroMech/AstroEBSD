function dataSet = h5saveRead(h5name, dataSetName,start,count,stride)
    % try to read data from h5 in parallel loops savely. we do this by waiting for 1000 iterations until the file is free again.
    % the file not being readable in parallel mode seems to be associated with a hdf5 library error.
    %
    %
    % L. Berners, 19.09.2025
    
    %% copying from matlab toolbox /matlab/imagesci/h5read.m
    switch(nargin)
        case 2
            start = [];
            count = [];
            stride = [];
        case 4
            stride = [];
        case 5    
            %
        otherwise
            error(message('MATLAB:imagesci:validate:wrongNumberOfInputs'));	   
    end

    % check whether we are in a parallel enviroentment
    isparfor = ~isempty(getCurrentTask());
    % disp(start)
    % disp(count)
    
    if ~isparfor
        dataSet= h5read(h5name,dataSetName,start,count,stride);
        return;
    end
    

    % try to read 1000 times, if no success we still throw an erro. but this way 
    % we might have a chance to read a file in a parallel loop.
    whileCount = 0;
    while 1
        try
            isparfor = 0;
            dataSet= h5read(h5name,dataSetName,start,count,stride);
        catch e
            if contains( e.message, "HDF5 Library error")
                isparfor = 1;
                pause(0.2)
            end
            warning('Dataset %s could not be read',dataSetName)
        end    
        
        if isparfor ==0
            break
        end

        if whileCount>1000
            break;
        end
        whileCoun = whileCount+1
    end
end