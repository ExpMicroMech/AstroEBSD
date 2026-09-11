function  InputUserOld = overwriteInputUser(InputUserOld, InputUserNew)
    fieldNames = fieldnames(InputUserNew);

    for fn= 1:size(fieldNames,1)
        InputUserOld.(fieldNames{fn}) = InputUserNew.(fieldNames{fn}) ;
    end
end