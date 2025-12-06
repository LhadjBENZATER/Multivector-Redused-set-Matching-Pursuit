function stop(~, ~)
    handles = guidata(gcf);
    handles.isRunning = false;

    % Close UDP connection if it exists
    if isfield(handles, 'udpObj') && ~isempty(handles.udpObj)
        clear handles.udpObj;
    end

    guidata(gcf, handles);
end
