function txt = DisplayCursor(Target, eventObj, vAx, mFourthVariable)
    if numel(vAx) == 1
        fourthVariable = mFourthVariable;
    else
        fourthVariable = mFourthVariable(:,Target.Parent == vAx);
    end
    % Customizes text of data tips, showing additional fourth variable
    pos = get(eventObj, 'Position');
    idx = get(eventObj, 'DataIndex');
    txt = {['X: ',num2str(pos(1))],...
           ['Y: ',num2str(pos(2))],...
           ['Z: ',num2str(pos(3))],...
           ['D: ',num2str(fourthVariable(idx))]};    
end