function trk = fix_GetPartOfTrack(trk,f0,f1)
% return subset of the input trx from frame f0 to f1

i0 = trk.off+(f0);
i1 = trk.off+(f1);
i0 = max(1,i0);
i1 = min(i1,trk.nframes);
idx = i0:i1;

nfrm = trk.nframes;
assert(numel(trk.timestamps)==nfrm,...
  'Invalid .timeseries field encountered in trx.');

fldsTS = Fix.FLDS_TIMESERIESDATA;
fldsUnk = setdiff(fieldnames(trk),Fix.FLDS_ALLKNOWN);
for f=fldsTS(:)',f=f{1}; %#ok<FXSET>
  if isfield(trk,f)
    trk.(f) = trk.(f)(idx);
  end
end
for f=fldsUnk(:)',f=f{1}; %#ok<FXSET>
  v = trk.(f);
  if isvector(v) && numel(v)==nfrm
    warningNoTrace('trk:unk','Unexpected timeseries field: %s',f);
    trk.(f) = trk.(f)(idx);
  else
    warningNoTrace('trk:unk','Unrecognized trx field: %s',f);
  end
end

trk.nframes = max(0,i1-i0+1);
trk.firstframe = max(f0,trk.firstframe);
trk.endframe = min(trk.endframe,f1);
trk.off = -trk.firstframe + 1;
assert(isequal(trk.nframes,trk.endframe-trk.firstframe+1,i1-i0+1));
