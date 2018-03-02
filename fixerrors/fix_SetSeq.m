function handles = fix_SetSeq(handles,seqi,isfirstframe)
% set the GUI state for displaying a particular sequence index
% splintered from fixerrorsgui 6/23/12 JAB

if exist('isfirstframe','var')==0
  isfirstframe = false;
end

handles.seqi = seqi;
seq = handles.seqs(seqi);
handles.f = seq.frames(1);
handles.nselect = 0;
handles.selected = [];
set(handles.errnumbertext,'string',sprintf('Error: %d/%d',seqi,length(handles.seqs)));
set(handles.seqframestext,'string',sprintf('Frames: %d:%d',seq.frames(1),seq.frames(end)));
set(handles.seqfliestext,'string',['Flies: [',num2str(seq.flies),']']);
set(handles.seqtypetext,'string',sprintf('Type: %s',seq.type));
set(handles.seqsusptext,'string',sprintf('Susp: %f',max(seq.suspiciousness)));

if ~isfirstframe
  seqTable = handles.seqTable;
  assert(seqTable.nRows==numel(handles.seqs));
  seqTable.setSelectedRows(seqi);
end

% AL20180302: unnec for JS and requires stats toolbox
% % set fly colors so that flies that are close have different colors
% x = nan(1,handles.nflies);
% y = nan(1,handles.nflies);
% f = round(mean([seq.frames(1),seq.frames(end)]));
% for fly = 1:handles.nflies,
%   if ~isalive(handles.trx(fly),f),
%     continue;
%   end
%   i = handles.trx(fly).off+(f);
%   x(fly) = handles.trx(fly).x(i);
%   y(fly) = handles.trx(fly).y(i);
% end

% D = squareform(pdist([x;y]'));
% handles.colors(seq.flies,:) = handles.colors0(handles.colororder(1:length(seq.flies)),:);
% isassigned = false(1,handles.nflies);
% isassigned(seq.flies) = true;
% D(:,seq.flies) = nan;
% for i = length(seq.flies)+1:handles.nflies,
%   [mind,fly] = min(min(D(isassigned,:),[],1));
%   if isnan(mind),
%     handles.colors(~isassigned,:) = handles.colors0(handles.colororder(i:end),:);
%     break;
%   end
%   handles.colors(fly,:) = handles.colors0(handles.colororder(i),:);
%   isassigned(fly) = true;
%   D(:,fly) = nan;
% end

if isfield(handles,'hpath'),
  for fly = 1:handles.nflies,
    if length( handles.hpath ) < fly
      fprintf( 1, 'error at fly %d: nflies %d; len hpath %d, len hcenter %d\n', fly, handles.nflies, length( handles.hpath ), length( handles.hcenter ) );
      break
    end
    clr = handles.colors(fly,:);
    safeset(handles.hpath(fly),'color',clr);
    safeset(handles.hpath(fly),'color',clr);
    safeset(handles.htailmarker(fly),'color',clr);
    safeset(handles.hellipse(fly),'color',clr);
    safeset(handles.hleft(fly),'color',clr);
    safeset(handles.hright(fly),'color',clr);
    safeset(handles.hhead(fly),'color',clr);
    safeset(handles.htail(fly),'color',clr);
    safeset(handles.hcenter(fly),'color',clr);
    safeset(handles.hwingl(fly),'color',clr);
    safeset(handles.hwingr(fly),'color',clr);
    safeset(handles.hwinglinel(fly),'color',clr);
    safeset(handles.hwingliner(fly),'color',clr);
  end
end

if ~isfirstframe
  fix_SetFrameNumber(handles);
  fix_PlotFrame(handles);
  fix_ZoomInOnSeq(handles);
end


function safeset(h,varargin)

if ishandle(h),
  set(h,varargin{:});
end
