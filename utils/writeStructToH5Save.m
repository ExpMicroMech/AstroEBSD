function writeStructToH5Save(struct, filename,path)

    % we use parfor both to check, whether run in parallel as well to wait
    % until data is read. assuming, that at some point the h5 must be fre
    % we should not wait longer than a certain time.

    isparfor = ~isempty(getCurrentTask());

    if ~isparfor
        writeStructToH5(struct, filename,path);
        return;
    end


    whileCount = 0;

    while 1
        try
            isparfor = 0;
            writeStructToH5(struct, filename,path);
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

        if whileCount>1000
            break;
        end
        whileCoun = whileCount+1
end
end