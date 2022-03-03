function UpdateCursorDataTip(fig, vAx, fourthVariable)
dcmObj = datacursormode(fig);
set(dcmObj, 'UpdateFcn', {@DisplayCursor, vAx, fourthVariable})
end

