function fix_Play(handles,hObject,speedfac)
% play through a sequence
% splintered from fixerrorsgui 6/23/12 JAB

if exist('speedfac','var')==0
  speedfac = 1.0;
end

handles.isplaying = true;
origstr = get(hObject,'string');
set(hObject,'string','Stop','backgroundcolor',[.5,0,0]);
guidata(hObject,handles);

switch get(hObject,'tag')
  case {'playstopbutton' 'playstopbuttonslow'}
    f0 = max(1,handles.seq.frames(1)-10);
    f1 = min(handles.nframes,handles.seq.frames(end)+10);
  otherwise
    f0 = handles.f;
    f1 = handles.nframes;
end     

minSPFuse = handles.MinSPF/speedfac;
tic;
for f = f0:f1  
  handles = guidata(hObject);
  if ~handles.isplaying
    break;
  end  
  
  handles.f = f;
  fix_SetFrameNumber(handles);
  fix_PlotFrame(handles);
  drawnow;
  handles = guidata(hObject);

  if handles.MaxFPS > 0
    tmp = toc;
    if tmp < minSPFuse
      pause(minSPFuse - tmp);
    end
  else
    drawnow;
  end
  tic;
end

if handles.isplaying
  handles.f = handles.seq.frames(1);
  fix_SetFrameNumber(handles);
  fix_PlotFrame(handles);  
end

handles.isplaying = false;
set(hObject,'string',origstr,'backgroundcolor',[0,.5,0]);
guidata(hObject,handles);
