function MyMsgBox(msgBoxMsg, msgBoxTitle, b_waitForOk)
if ~exist('b_waitForOk', 'var')
    b_waitForOk = false;
end
sMsgBox.Interpreter = 'tex';
sMsgBox.WindowStyle = 'modal';
msgBoxIcon = 'warn';
if b_waitForOk
    answerStr = questdlg(msgBoxMsg, msgBoxTitle);
    if strcmp(answerStr,'Yes')
        answer = true;
    else
        answer = false;
    end
else
    msgBoxMsg = strcat('\fontsize{12} ',msgBoxMsg);
    msgbox(msgBoxMsg, msgBoxTitle, msgBoxIcon, sMsgBox)
    answer = true;
end
assert(answer);
end