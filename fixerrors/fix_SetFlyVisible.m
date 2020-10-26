function fix_SetFlyVisible(handles,fly,v)
% makes fly body visible or invisible
% splintered from fixerrorsgui 6/21/12 JAB
% Modified by Allen Lee for DTFE in 2017

if isdummytrk(handles.trx(fly))
  return;
end

set(handles.hellipse(fly),'visible',v);
set(handles.hcenter(fly),'visible',v);
set(handles.hleft(fly),'visible',v);
set(handles.hright(fly),'visible',v);
set(handles.hhead(fly),'visible',v);
set(handles.htail(fly),'visible',v);
set(handles.htailmarker(fly),'visible',v);
set(handles.hpath(fly),'visible','off');
set(handles.hwingl(fly),'visible',v);
set(handles.hwingr(fly),'visible',v);
set(handles.hwinglinel(fly),'visible',v);
set(handles.hwingliner(fly),'visible',v);
