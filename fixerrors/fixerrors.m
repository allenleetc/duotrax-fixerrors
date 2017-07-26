function fixerrors

%% path
mpath = mfilename('fullpath');
mpath = fileparts(mpath);
addpath(mpath);
addpath(fullfile(mpath,'JavaTableWrapper'));
setuppath;

%% Movie
MOVIEEXTS = {'*.avi','*.fmf','*.sbfmf','*.ufmf'}';
HELPMSG = 'Choose movie file for which to fix tracking errors';
moviename = Fix.cfgGetProp('moviename');
if isempty(moviename)
  moviename = pwd;
end
[movienameS,moviepath] = uigetfilehelp(MOVIEEXTS,'Choose movie file',...
  moviename,'helpmsg',HELPMSG);
if isequal(movienameS,0)
  return;
end
moviename = fullfile(moviepath,movienameS);
Fix.cfgSetProp('moviename',moviename);

%% Look for restart file
savedProgFiles = Fix.findSavedProgFiles(moviename);
tfRestart = false;
restartFile = [];
restartContents = [];
if ~isempty(savedProgFiles)
  liststr = [{'<Don''t restart, start fresh>'};savedProgFiles(:)];
  [sel,ok] = listdlg('ListString',liststr,...
    'SelectionMode','single',...
    'Name','Restart saved progress',...
    'PromptString','Saved progress found for this movie. Restart?',...
    'OKString','Restart with selected file');
  if ok==0
    return;
  end
  tfRestart = sel>1;
  if tfRestart
    restartFile = liststr{sel};
    restartContents = load(restartFile);
    fprintf(1,'Loaded restart file %s...\n',restartFile);
    if ~strcmp(restartContents.moviename,moviename)
      error('fixerrors:restart',...
        'Moviename in restart file (%s) does not match selected movie (%s).',...
        restartContents.moviename,moviename);
    end
    
    % Confirm restart with stored trxfile
    qstr = sprintf('Restart file contents:\nmovie: %s\ntrx: %s\n.',...
      restartContents.moviename,restartContents.matname);
    resp = questdlg(qstr,'Confirm restart','OK, Restart','Cancel','OK, Restart');
    if isempty(resp)
      resp = 'Cancel';
    end
    switch resp
      case 'OK, Restart'
        % none, tfRestart, restartFile, restartContents all set.
      case 'Cancel'
        return;
      otherwise
        assert(false);
    end
  end
end

%% Set: seqs, trx, annname, params, trxname, trxnameS, undolist
if tfRestart
  seqs = restartContents.seqs;
  trx = restartContents.trx;
  annname = restartContents.annname;
  params = restartContents.params;
  trxname = restartContents.matname;
  undolist = restartContents.undolist;
  
  [~,trxnameS,ext] = fileparts(trxname);
  trxnameS = [trxnameS ext];
  
  % legacy
  assert(~isfield(trx,'f2i'));
  assert(isfield(trx,'off'));  
else
    
  %% Trx
  HELPMSG = sprintf('Choose trx file containing trajectories for movie %s.',...
    moviename);
  trxname = fullfile(moviepath,'registered_trx.mat');
  [trxnameS,trxpath] = uigetfilehelp({'*.mat'},'Choose trx file',trxname,...
    'helpmsg',HELPMSG);
  if isequal(trxnameS,0)
    return;
  end
  trxname = fullfile(trxpath,trxnameS);
  annname = fullfile(trxpath,'READANN_DUMMY');

  [seqs,trx,params] = suspicious_sequences(trxname,annname,...
    'minerrjumpfrac',nan,'minorientchange',nan,...
    'maxmajorfrac',nan,'minwalkvel',nan,...
    'matcherrclose',nan,'minanglediff',nan);
  undolist = {};
  if isempty(seqs)
    msgbox('No suspicious sequences found.');
    return;
  end

%   [max_jump,maxmajor,meanmajor] = read_ann(annname,...
%     'max_jump','maxmajor','meanmajor');
%   meanmajor = meanmajor * 4;
%   maxmajor = maxmajor * 4;
  
%   px2mm = trx(1).pxpermm; % pixels per millimeter
%   Fix.cfgSetProp('px2mm',px2mm);  
%   max_jump = max_jump / px2mm;
%   maxmajor = maxmajor / px2mm;
%   meanmajor = meanmajor / px2mm;
%   mm2body = meanmajor; % millimeters per body-length
  
%   % set default values
%   if ~exist('minerrjump','var')
%     minerrjump = .2*max_jump;
%   end
%   if ~exist('minorientchange','var')
%     minorientchange = nan;
%   end
%   if ~exist('largemajor','var')
%     largemajor = meanmajor + 2/3*(maxmajor-meanmajor);
%   end
%   if ~exist('minanglediff','var')
%     minanglediff = nan;
%   end
%   if ~exist('minwalkvel','var')
%     minwalkvel = 1 / 4;
%   end
%   if ~exist('matcherrclose','var')
%     matcherrclose = 10/4^2;
%   end
%   tmp = [minerrjump/mm2body, ...
%     minorientchange, ...
%     largemajor, ...
%     minanglediff, ...
%     minwalkvel/mm2body/trx(1).fps, ...
%     matcherrclose/mm2body/mm2body];
%   defaultv = cell(size(tmp));
%   for i = 1:length(tmp)
%     defaultv{i} = num2str( tmp(i) );
%   end
%   
%   shortdescr = cell(1,6);
%   descr = cell(1,6);
%   relev = cell(1,6);
%   shortdescr{1} = 'Value 1: Minimum suspicious prediction error (body-lengths)';
%   descr{1} = ['All sequences in which the error between the constant-velocity ',...
%     'prediction and the fly''s measured position is greater than Value 1 ',...
%     'will be flagged.'];
%   relev{1} = sprintf('Max. jump error: %.1f mm; Mean body length: %.1f mm',max_jump, mm2body);
%   shortdescr{2} = 'Value 2: Minimum suspicious orientation change (deg)';
%   descr{2} = ['All sequences in which the change in orientation is greater ',...
%     'than Value 2 will be flagged.'];
%   relev{2} = '';
%   shortdescr{3} = 'Value 3: Minimum suspiciously large major axis (mm)';
%   descr{3} = ['All sequences in which the major axis length (i.e., the body length) ',...
%     'is greater than Value 3 will be flagged.'];
%   relev{3} = sprintf('Mean major axis length: %.2f mm; max. major axis length: %.2f mm',meanmajor,maxmajor);
%   shortdescr{4} = 'Value 4a: Minimum suspicious orientation-direction mismatch (deg): ';
%   descr{4} = '';
%   relev{4} = '';
%   shortdescr{5} = 'Value 4b: Minimum walking speed (body-lengths/sec)';
%   descr{5} = ['All sequences in which the fly is walking (has speed greater than ',...
%     'Value 4b) and its body orientation differs from its movement direction ',...
%     'by more than Value 4a value will be flagged.'];
%   relev{5} = sprintf( 'Mean body length: %.1f mm; frames/sec: %.f', mm2body, trx(1).fps );
%   shortdescr{6} = 'Value 5: Maximum ambiguous-identity error (body-lengths^2)';
%   descr{6} = ['All sequences in which the error added by swapping ',...
%     'the identities of two fly tracks is less than Value 5 will be flagged.'];
%   relev{6} = sprintf( 'Mean body length: %.1f mm', mm2body );
%   assert( length( defaultv ) == length( shortdescr ) )
%   prompts = cell(size(shortdescr));
%   for i = 1:length(shortdescr)
%     prompts{i} = sprintf('**%s**: %s',shortdescr{i},descr{i});
%     if ~isempty(relev{i})
%       prompts{i} = [prompts{i},sprintf(' [Relevant quantities: %s]',relev{i})];
%     end
%   end
%   title1 = 'Suspiciousness Parameters';
%   tmp = inputdlg(prompts,title1,1,defaultv,'on');
%   
%   if isempty(tmp)
%     return;
%   end
  
%   minerrjump = str2double(tmp{1})*mm2body; % convert back to mm
%   minorientchange = str2double(tmp{2});
%   largemajor = str2double(tmp{3});
%   minanglediff = str2double(tmp{4});
%   minwalkvel = str2double(tmp{5})*mm2body*trx(1).fps;
%   matcherrclose = str2double(tmp{6})*mm2body*mm2body;
  
%   try
%     save('-append',savedsettingsfile,...
%       'minerrjump','minorientchange','largemajor','minanglediff',...
%       'minwalkvel','matcherrclose');
%   catch ME
%     fprintf('Could not save to settings file %s -- not a big deal\n',savedsettingsfile);
%     getReport(ME)
%   end
  
  %% convert to the units expected by suspicious_sequences
%   minerrjumpfrac = minerrjump / max_jump;
%   minorientchange = minorientchange*pi/180;
%   maxmajorfrac = (largemajor - meanmajor)/(maxmajor - meanmajor);
%   minwalkvel = minwalkvel*px2mm;
%   matcherrclose = matcherrclose*px2mm^2;
%   minanglediff = minanglediff*pi/180;
%   [seqs,trx0,params] = suspicious_sequences(trxname,annname,...
%     'minerrjumpfrac',minerrjumpfrac,'minorientchange',minorientchange,...
%     'maxmajorfrac',maxmajorfrac,'minwalkvel',minwalkvel,...
%     'matcherrclose',matcherrclose,'minanglediff',minanglediff);
end

%% Call fixerrorsgui

saveProgFile = Fix.createSavedProgFilename(movienameS,trxnameS);
saveProgFile = fullfile(moviepath,saveProgFile);
fprintf('Movie: %s\n',moviename);
fprintf('Trx: %s\n',trxname);
fprintf('Temporary progress file will be created at: %s\n',saveProgFile);

if ~isfield(seqs,'status')
  [seqs.status] = deal(SeqStatus.UNKNOWN);
end

readfrm = struct();
[readfrm.readframe,readfrm.nframes,readfrm.fid] = get_readframe_fcn(moviename);
fixerrorsgui(seqs,moviename,trx,annname,params,trxname,saveProgFile,...
  readfrm,undolist);
