function writeStructToH5(struct, filename,path)

% Partially generate using RWTH gpt (OPENAIGPT 40 mini)

    fields = fieldnames(struct);
    for i = 1:length(fields)
        field = fields{i};
        value = struct.(field);
        
        % Create a new path for this field
        newPath = sprintf("/%s/%s",path,field);
        
        if isstruct(value) % If it's a nested struct, recurse
            writeStructToH5(value, filename,newPath);
        else 
            % Determine the size of the value to create the dataset appropriately
            if isnumeric(value) || islogical(value)
               try 
                   h5create(filename, newPath, size(value));
               catch
                   warning('DAtaset allready exists')
               end
                 try
                        h5write(filename, newPath, value);
                    catch
                        warning(sprintf("Dataset %s Value could not be written",newPath))
                    end
                %%
            elseif ischar(value) || isstring(value) % Handle strings separately
                % h5create(filename, newPath, size(value));
                try 
                   h5create(filename, newPath, size(string(value)),'Datatype','string');
               catch
                   warning('DAtaset allready exists')
               end
                     try
                        h5write(filename, newPath, string(value));
                    catch
                        warning(sprintf("Dataset %s Value could not be written",newPath))
                    end
            elseif iscell(value)
                %%
                for cellNum = 1:size(value(:),1)
                    % outsurce this to writeCellToH5 in the future.
                    try 
                       h5create(filename, sprintf("%s/%i",newPath,cellNum), size(string(value{cellNum})),'Datatype','string');
                   catch
                       warning('DAtaset allready exists')
                    end
                    try
                        h5write(filename, sprintf("%s/%i",newPath,cellNum), string(value{cellNum}));
                    catch
                        warning(sprintf("Dataset %s Value could not be written",newPath))
                        % warning("Dataset Value could not be written")
                    end
                end
                
            % end
            else
                warning('Unsupported data type for field: %s', field);
            end
        end
    end
end