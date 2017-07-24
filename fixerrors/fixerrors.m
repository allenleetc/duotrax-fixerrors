function fixerrors

% prompts user for mat, annotation, and movie files, parameters
% for computing suspicious frames, then computes suspicious frames, then
% brings up the fixerrors gui

%% path
mpath = mfilename('fullpath');
mpath = fileparts(mpath);
addpath(mpath);
addpath(fullfile(mpath,'JavaTableWrapper'));

%% set all defaults

moviename = '';
moviepath = '';
setuppath;

%% read last settings
pathtofixerrors = which('fixerrors');
savedsettingsfile = strrep(pathtofixerrors,'fixerrors.m','.fixerrorsrc.mat');
if exist(savedsettingsfile,'file')
  load(savedsettingsfile);
end

%% Movie
fprintf('Choose a movie to fix errors in\n\n');
movieexts = {'*.avi','*.fmf','*.sbfmf','*.ufmf'}';
helpmsg = 'Choose movie file for which to fix tracking errors';
[moviename,moviepath] = uigetfilehelp(movieexts,'Choose movie file',moviename,'helpmsg',helpmsg);
if isnumeric(moviename) && moviename==0
  return;
end
[movietag,movieext] = splitext(moviename);

%% Trx
helpmsg = sprintf('Choose the trx file containing the trajectories corresponding to movie %s.',[moviepath,moviename]);
matname = [moviepath,strrep(moviename,movieext,'.mat')];
[matname,matpath] = uigetfilehelp({'*.mat'},'Choose trx file',matname,'helpmsg',helpmsg);
if isnumeric(matname) && matname == 0
  return;
end

%% Ann
annname = fullfile(matpath,'READANN_DUMMY');
% annname = [matpath,'roidata.mat'];
% helpmsg = {};
% helpmsg{1} = 'Choose the ROI data file corresponding to:';
% helpmsg{2} = sprintf('Movie: %s; and',[moviepath,moviename]);
% helpmsg{3} = sprintf('Trajectory mat file: %s',[matpath,matname]);
% [annname,annpath] = uigetfilehelp({'*.mat'},'Choose roidata file',annname,'helpmsg',helpmsg);
% if isnumeric(annname) && annname == 0
%   return;
% end

%% RC file
tmpfilename = fullfile(matpath,sprintf('fixerrorsTmpProg_%s_%s.mat',...
  splitext(moviename),splitext(matname)));
moviename = [moviepath,moviename];
matname = [matpath,matname];
%annname = [annpath,annname];
try
  if exist('savedsettingsfile','file')
    save('-append',savedsettingsfile,'moviename','moviepath');
  else
    save(savedsettingsfile,'moviename','moviepath');
  end
catch ME
  fprintf('Could not save to rc file %s -- not a big deal\n',savedsettingsfile);
  getReport(ME)
end

%% 
[readframe,nframes,fid] = get_readframe_fcn(moviename);
readframe_fcn.readframe = readframe;
readframe_fcn.nframes = nframes;
readframe_fcn.fid = fid;

%% Restart
tag = movietag;
DORESTART = false;

if exist(tmpfilename,'file')
  tmp2 = load(tmpfilename,'matname','moviename');
  if isfield(tmp2,'moviename')
    oldmoviename = tmp2.moviename;
  else
    oldmoviename = 'unknown movie';
  end
  if isfield(tmp2,'matname')
    oldmatname = tmp2.matname;
  else
    oldmatname = 'unknown trx file';
  end
  prompt = {};
  prompt{1} = sprintf('A restart file saved by fixerrors was found with tag %s ',tag);
  prompt{2} = sprintf('Original movie: %s, selected movie %s. ',oldmoviename,moviename);
  prompt{3} = sprintf('Original trx file: %s, selected trx file %s. ',oldmatname,matname);
  prompt{4} = 'It is only recommended that you load these partial results if you are certain the trx files match. ';
  prompt{5} = 'Would you like to load the saved results and restart? ';
  
  button = questdlg(cell2mat(prompt),'Restart?');
  if strcmpi(button,'yes')
    DORESTART = true;
  elseif strcmpi(button,'cancel')
    return
  end
end


%% convert to px, seconds
trx = load(matname);
trx = trx.trx;
% [matpathtmp,matnametmp] = split_path_and_filename(matname);
% [convertsucceeded,convertmatname,trx] = convert_units_f('matname',matnametmp,...
%    'matpath',matpathtmp,'moviename',moviename,'readframe_fcn',readframe_fcn);
% if ~convertsucceeded,
%   return;
% end
% convertmatname = matname;
%[trx,matname,succeeded] = load_tracks(convertmatname,moviename);

%% set parameters for detecting suspicious frames

if ~DORESTART
  
  [max_jump,maxmajor,meanmajor,arena_radius, arena_center_x, arena_center_y, do_set_circular_arena] = ...
    read_ann(annname,'max_jump','maxmajor','meanmajor','arena_radius', 'arena_center_x', 'arena_center_y', 'do_set_circular_arena');
  meanmajor = meanmajor * 4; % why 4?
  maxmajor = maxmajor * 4;
  
  circular_arena = struct( ...
    'do_set_circular_arena', do_set_circular_arena, ...
    'arena_radius', arena_radius, ...
    'arena_center_x', arena_center_x, ...
    'arena_center_y', arena_center_y );
  
  px2mm = trx(1).pxpermm; % pixels per millimeter
  try
    save('-append',savedsettingsfile,'px2mm');
  catch ME
    fprintf('Could not save to settings file %s -- not a big deal\n',savedsettingsfile);
    getReport(ME)
  end
  
  max_jump = max_jump / px2mm;
  maxmajor = maxmajor / px2mm;
  meanmajor = meanmajor / px2mm;
  mm2body = meanmajor; % millimeters per body-length
  
  % set default values
  if ~exist('minerrjump','var')
    minerrjump = .2*max_jump;
  end
  if ~exist('minorientchange','var')
    minorientchange = nan;
  end
  if ~exist('largemajor','var')
    largemajor = meanmajor + 2/3*(maxmajor-meanmajor);
  end
  if ~exist('minanglediff','var')
    minanglediff = nan;
  end
  if ~exist('minwalkvel','var')
    minwalkvel = 1 / 4;
  end
  if ~exist('matcherrclose','var')
    matcherrclose = 10/4^2;
  end
  tmp = [minerrjump/mm2body, ...
    minorientchange, ...
    largemajor, ...
    minanglediff, ...
    minwalkvel/mm2body/trx(1).fps, ...
    matcherrclose/mm2body/mm2body];
  defaultv = cell(size(tmp));
  for i = 1:length(tmp)
    defaultv{i} = num2str( tmp(i) );
  end
  
  shortdescr = cell(1,6);
  descr = cell(1,6);
  relev = cell(1,6);
  shortdescr{1} = 'Value 1: Minimum suspicious prediction error (body-lengths)';
  descr{1} = ['All sequences in which the error between the constant-velocity ',...
    'prediction and the fly''s measured position is greater than Value 1 ',...
    'will be flagged.'];
  relev{1} = sprintf('Max. jump error: %.1f mm; Mean body length: %.1f mm',max_jump, mm2body);
  shortdescr{2} = 'Value 2: Minimum suspicious orientation change (deg)';
  descr{2} = ['All sequences in which the change in orientation is greater ',...
    'than Value 2 will be flagged.'];
  relev{2} = '';
  shortdescr{3} = 'Value 3: Minimum suspiciously large major axis (mm)';
  descr{3} = ['All sequences in which the major axis length (i.e., the body length) ',...
    'is greater than Value 3 will be flagged.'];
  relev{3} = sprintf('Mean major axis length: %.2f mm; max. major axis length: %.2f mm',meanmajor,maxmajor);
  shortdescr{4} = 'Value 4a: Minimum suspicious orientation-direction mismatch (deg): ';
  descr{4} = '';
  relev{4} = '';
  shortdescr{5} = 'Value 4b: Minimum walking speed (body-lengths/sec)';
  descr{5} = ['All sequences in which the fly is walking (has speed greater than ',...
    'Value 4b) and its body orientation differs from its movement direction ',...
    'by more than Value 4a value will be flagged.'];
  relev{5} = sprintf( 'Mean body length: %.1f mm; frames/sec: %.f', mm2body, trx(1).fps );
  shortdescr{6} = 'Value 5: Maximum ambiguous-identity error (body-lengths^2)';
  descr{6} = ['All sequences in which the error added by swapping ',...
    'the identities of two fly tracks is less than Value 5 will be flagged.'];
  relev{6} = sprintf( 'Mean body length: %.1f mm', mm2body );
  assert( length( defaultv ) == length( shortdescr ) )
  prompts = cell(size(shortdescr));
  for i = 1:length(shortdescr)
    prompts{i} = sprintf('**%s**: %s',shortdescr{i},descr{i});
    if ~isempty(relev{i})
      prompts{i} = [prompts{i},sprintf(' [Relevant quantities: %s]',relev{i})];
    end
  end
  title1 = 'Suspiciousness Parameters';
  tmp = inputdlg(prompts,title1,1,defaultv,'on');
  
  if isempty(tmp)
    return;
  end
  
  minerrjump = str2double(tmp{1})*mm2body; % convert back to mm
  minorientchange = str2double(tmp{2});
  largemajor = str2double(tmp{3});
  minanglediff = str2double(tmp{4});
  minwalkvel = str2double(tmp{5})*mm2body*trx(1).fps;
  matcherrclose = str2double(tmp{6})*mm2body*mm2body;
  
  try
    save('-append',savedsettingsfile,...
      'minerrjump','minorientchange','largemajor','minanglediff',...
      'minwalkvel','matcherrclose');
  catch ME
    fprintf('Could not save to settings file %s -- not a big deal\n',savedsettingsfile);
    getReport(ME)
  end
  
  %% convert to the units expected by suspicious_sequences
  minerrjumpfrac = minerrjump / max_jump;
  minorientchange = minorientchange*pi/180;
  maxmajorfrac = (largemajor - meanmajor)/(maxmajor - meanmajor);
  minwalkvel = minwalkvel*px2mm;
  matcherrclose = matcherrclose*px2mm^2;
  minanglediff = minanglediff*pi/180;
  [seqs,trx0,params] = suspicious_sequences(matname,annname,...
    'minerrjumpfrac',minerrjumpfrac,'minorientchange',minorientchange,...
    'maxmajorfrac',maxmajorfrac,'minwalkvel',minwalkvel,...
    'matcherrclose',matcherrclose,'minanglediff',minanglediff);
  
end

%% call the fixerrors gui

fprintf('Movie: %s\n',moviename);
fprintf('Mat: %s\n',matname);
%fprintf('Annname: %s\n',annname);
fprintf('Temporary progress file will be created at: %s\n',tmpfilename);

if DORESTART
  realmatname = matname;
  load(tmpfilename);
  if isfield(trx,'f2i')
    trx = rmfield(trx,'f2i');
  end
  if ~isfield(trx,'off')
    for i = 1:length(trx)
      trx(i).off = -trx(i).firstframe + 1;
    end
  end
  matname = realmatname;
  trx0 = trx;
  clear trx;
  
  if ~exist( 'circular_arena', 'var' )
    circular_arena = struct( 'do_set_circular_arena', 0 );
  end
else
  undolist = {};
end

if ~isfield(seqs,'status')
  [seqs.status] = deal(SeqStatus.UNKNOWN);
end

trx = fixerrorsgui(seqs,moviename,trx0,annname,params,matname,tmpfilename,...
  readframe_fcn,circular_arena,undolist);
