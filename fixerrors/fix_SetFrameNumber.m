function fix_SetFrameNumber(handles,hObject)
% show the current frame
% splintered from fixerrorsgui 6/21/12 JAB

if nargin < 2,
  hObject = -1;
end
% if ~isfield( handles, 'seq' ) % maybe during initialization
%    fprintf( 1, 'ignoring uninitialized handles during frame slider callback\n' );
%    return 
% end

if hObject ~= handles.frameslider,
  set(handles.frameslider,'Value',handles.f);
end
if hObject ~= handles.frameedit,
  set(handles.frameedit,'string',num2str(handles.f));
end
seq = handles.seqs(handles.seqi);
if handles.f < seq.frames(1),
  set(handles.frameofseqtext,'string','Before Sequence','backgroundcolor',[1,0,0],...
    'foregroundcolor',[1,1,1]);
elseif handles.f > seq.frames(end),
  set(handles.frameofseqtext,'string','After Sequence','backgroundcolor',[1,0,0],...
    'foregroundcolor',[1,1,1]);
elseif handles.f == seq.frames(1),
  set(handles.frameofseqtext,'string','Frame of Seq: 1','backgroundcolor',[0,0,1],...
    'foregroundcolor',[1,1,1]);
elseif handles.f == seq.frames(end),
  set(handles.frameofseqtext,'string',...
    sprintf('Frame of Seq: %d',handles.f-seq.frames(1)+1),...
    'backgroundcolor',[1,1,0]/2,'foregroundcolor',[1,1,1]);
else
  set(handles.frameofseqtext,'string',...
    sprintf('Frame of Seq: %d',handles.f-seq.frames(1)+1),...
    'backgroundcolor',[.7,.7,.7],'foregroundcolor',[0,0,0]);
end
i = find(seq.frames == handles.f);
if isempty(i),
  set(handles.seqsusptext,'string','Susp: --');
else
  set(handles.seqsusptext,'string',sprintf('Susp: %f',seq.suspiciousness(i)));
end
% if handles.frmIsSeq(handles.f)
%   handles.txIsTouch.Visible = 'on';
% else
handles.txIsTouch.Visible = 'off';
% end
if ~isempty(handles.motionobj),
  if ~isalive(handles.trx(handles,motionobj{2}),handles.f),
    handles.motionobj = [];
  end
end
