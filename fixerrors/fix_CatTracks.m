function trk1 = fix_CatTracks(trk1,trk2)
% concatenates trk2 onto trk1
% does not copy all fields -- convert_units must be re-run on the output track
% splintered from fixerrorsgui 6/21/12 JAB

flds = fieldnames(trk1);
assert(isequal(flds,fieldnames(trk2)));

tf = strncmp(flds,'susp',4);
fldsSusp = flds(tf);

n = trk2.nframes;

fldsTS = Fix.FLDS_TIMESERIESDATA;
fldsMatch = Fix.FLDS_MATCH;
fldsKnown = [Fix.FLDS_ALLKNOWN; fldsSusp];
fldsUnk = setdiff(flds,fldsKnown);
for f=fldsTS(:)',f=f{1}; %#ok<FXSET>
  if isfield(trk1,f)
    trk1.(f)(end+1:end+n) = trk2.(f);
  end
end
for f=fldsMatch(:)',f=f{1}; %#ok<FXSET>
  if ~isequaln(trk1.(f),trk2.(f))
    warningNoTrace('trk:fld','Expected trx fields ''%s'' to match.',f);
  end
end
for f=fldsUnk(:)',f=f{1}; %#ok<FXSET>
  v1 = trk1.(f);
  v2 = trk2.(f);
  if isvector(v1) && numel(v1)==trk1.nframes && ...
     isvector(v2) && numel(v2)==n
    warningNoTrace('trk:unk','Unexpected timeseries field: %s',f);
    trk1.(f)(end+1:end+n) = trk2.(f);
  else
    warningNoTrace('trk:unk','Unknown trk field: %s',f);
  end
end

trk1.nframes = trk1.nframes + n;
trk1.endframe = trk1.endframe+n;
