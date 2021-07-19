function [] = MyMsgBox(msgBoxMsg, msgBoxTitle, b_waitForOk)
if ~exist('b_waitForOk', 'var')
    b_waitForOk = false;
end
sMsgBox.Interpreter = 'tex';
sMsgBox.WindowStyle = 'modal';
msgBoxMsg = strcat('\fontsize{12} ',msgBoxMsg);
msgBoxIcon = 'warn';
if b_waitForOk
    uiwait(msgbox(msgBoxMsg, msgBoxTitle, msgBoxIcon, sMsgBox));
else
    msgbox(msgBoxMsg, msgBoxTitle, msgBoxIcon, sMsgBox)
end
end