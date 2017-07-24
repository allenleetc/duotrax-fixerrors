function fix_Play(handles,hObject,varargin)
% play through a sequence
% splintered from fixerrorsgui 6/23/12 JAB

[playdir,speedfac,speedfacseq] = myparse(varargin,...
  'playdir',1,...
  'speedfac',1.0,...
  'speedfacseq',1.0);

maxnonseq = isinf(speedfac);

handles.isplaying = true;
origstr = get(hObject,'string');
set(hObject,'string','Stop','backgroundcolor',[.5,0,0]);
guidata(hObject,handles);

switch get(hObject,'tag')
  case {'playstopbutton' 'playstopbuttonslow'}
    f0 = max(1,handles.seqs(handles.seqi).frames(1)-10);
    f1 = min(handles.nframes,handles.seqs(handles.seqi).frames(end)+10);
    switch playdir 
      case 1, frmsPlay=f0:f1;
      case -1, frmsPlay=f1:-1:f0;
    end
  otherwise
    switch playdir
      case 1, frmsPlay = handles.f:handles.nframes;
      case -1, frmsPlay = handles.f:-1:1;
    end
end     

minSPF = handles.MinSPF/speedfac;
minSPFseq = handles.MinSPF/speedfacseq;
tic;
for f = frmsPlay
  handles = guidata(hObject);
  if ~handles.isplaying
    break;
  end  
  
  handles.f = f;
  fix_SetFrameNumber(handles);
  fix_PlotFrame(handles);
  guidata(hObject,handles);
  drawnow;
  %handles = guidata(hObject);

  if handles.frmIsSeq(f)
    dtFrm = toc;
    if dtFrm < minSPFseq
      pause(minSPFseq - dtFrm);
    end
  else
    if ~maxnonseq
      dtFrm = toc;
      if dtFrm < minSPF
        pause(minSPF - dtFrm);
      end
    end
  end
  tic;
end

if handles.isplaying
  handles.f = handles.seqs(handles.seqi).frames(1);
  fix_SetFrameNumber(handles);
  fix_PlotFrame(handles);  
end

handles.isplaying = false;
set(hObject,'string',origstr,'backgroundcolor',[0,.5,0]);
guidata(hObject,handles);
