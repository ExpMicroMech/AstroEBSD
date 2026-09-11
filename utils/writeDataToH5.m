function writeDataToH5(filename,dataSetName,Data,Datatype)
    isparfor = ~isempty(getCurrentTask());
    % we use parfor both to check, whether run in parallel as well to wait
    % until data is read. assuming, that at some point the h5 must be fre
    % we should not wait longer than a certain time.
    whileCount = 0;
    while 1
        if nargin == 4
            try
                 h5create(filename,dataSetName,size(Data),'Datatype',Datatype)
            catch
                warning('Dataset %s allready exists',dataSetName)
            end
        else
            try
                 h5create(filename,dataSetName,size(Data))
              catch
                warning('Dataset %s allready exists',dataSetName)
            end
        end
        try
            isparfor = 0;
            h5write(filename, dataSetName,Data)
        catch e
            if contains( e.message, "HDF5 Library error")
                isparfor = 1;
                pause(0.2)
            end
            warning('Dataset %s could not be written',dataSetName)
        end    
        if isparfor ==0
            break
        end
        whileCount=whileCount+1
        if whileCount>1000
            break;
        end
    end
end