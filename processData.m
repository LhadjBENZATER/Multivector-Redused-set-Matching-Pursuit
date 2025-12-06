function processData(udpObj, ~)
    handles = guidata(gcf);
    if ~handles.isRunning
        return;
    end

    % Read the received data
    dataPacket = read(udpObj, udpObj.NumBytesAvailable, "double"); 

    if ~isempty(dataPacket)
        % Process data
        x = reshape(dataPacket, [], 1); % Reshape if needed
        x = fct_syn(x); % Apply your existing processing functions
        x = fct_Phase(x);

        % Store the data for further processing in the main loop
        handles.x = x; 
        guidata(gcf, handles);
    end
end
