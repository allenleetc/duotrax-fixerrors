function [tfdone,nexttype] = fix_SetErrorTypes(handles)
% * Find remaining uncorrected seqs in handles.seqs
% * Set nexterrortype menu based on remaining seqs
% * Set pbCorrect string based on done-ness
% 
% tfdone: scalar logical. if true, there are no remaining seqs
% nexttype: char. Next type to consider. Indeterminate if tfdone==true.

seqs = handles.seqs(:);
tf = [seqs.status]'==SeqStatus.UNKNOWN;
szassert(tf,size(seqs));
seqsRemain = seqs(tf);
typesRemain = unique({seqsRemain.type}','stable');
prettytypesRemain = cellfun(@SeqTypes.type2pretty,typesRemain,'uni',0);

pumNET = handles.nexterrortypemenu;
pbCorrect = handles.correctbutton;
selOrig = getListControlSelection(pumNET);

if isempty(prettytypesRemain)
  set(pumNET,'String','No more errors','Value',1,'UserData',[]);
  pbCorrect.Enable = 'off';
  %set(pbCorrect,'String','Finish');
  tfdone = true;
  nexttype = [];
else
  i = find(strcmpi(selOrig,prettytypesRemain));
  if ~isempty(i)
    val = i;
  else
    val = min(pumNET.Value,length(prettytypesRemain));
  end
  set(pumNET,'Value',val,'String',prettytypesRemain,'UserData',typesRemain);
  pbCorrect.Enable = 'on';
  %pbCorrect.String = 'Correct';
  tfdone = false;
  nexttype = typesRemain{val};
end
