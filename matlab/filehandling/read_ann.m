% [param1,param2,...] = read_ann(filename,'param1','param2',...)
%
% params:
%
% version, bg_type, n_bg_std_thresh, n_bg_std_thresh_low, bg_std_min,
% bg_std_max n_bg_frames, min_nonarena, max_nonarena, arena_center_x,
% arena_center_y, arena_radius, do_set_circular_arena, bg_algorithm,
% background_median, bg_norm_type, background_mad, hfnorm, bg_norm_type,
% hm_cutoff, hm_boost, hm_order, maxarea, maxmajor, maxminor, maxecc,
% minarea, minmajor, minminor, minecc, meanarea, meanmajor, meanminor,
% meanecc, nframes_size, nstd_shape, max_jump, ang_dist_wt, center_dampen,
% angle_dampen, minbackthresh, maxpenaltymerge, maxareadelete,
% do_fix_split, splitdetection_length, splitdetection_cost, do_fix_merged,
% mergeddetection_length, mergeddetection_distance, do_fix_spurious,
% spuriousdetection_length, do_fix_lost, lostdetection_length, movie_name,
% start_frame, data_format, velocity_angle_weight,
% max_velocity_angle_weight, 
% fracframesisback, expbgfgmodel_filename, use_expbgfgmodel, 
% expbgfgmodel_llr_thresh, min_frac_frames_isback,
% background_center, background_dev, movie_height, movie_width

% Modified by Allen Lee for DTFE in 2017

function varargout = read_ann(filename,varargin)

[p,f,e] = fileparts(filename);
filenameS = [f e];
if strcmp(filenameS,'READANN_DUMMY')
  roifname = fullfile(p,'roidata.mat');
  tdfname = fullfile(p,'trackingdata.mat');
  bgfname = fullfile(p,'bgdata.mat');
  %fprintf(1,'Fake read_ann from %s/%s/%s.\n',roifname,tdfname,bgfname);  
  roi = load(roifname);
  td = load(tdfname);
  bg = load(bgfname);
  tdparams = td.params;
  
  s = struct();
  s.center_dampen = tdparams.err_dampen_pos;  
  s.angle_dampen = tdparams.err_dampen_theta;
  s.max_jump = nan;
  s.maxmajor = nan;
  s.meanmajor = nan;
  s.velocity_angle_weight = nan; % for ID
  s.ang_dist_wt = nan; % for ID
  s.arena_radius = roi.radii;
  s.arena_center_x = roi.centerx;
  s.arena_center_y = roi.centery;
  s.do_set_circular_arena = 1;
  s.bg_algorithm = 'median';
  s.background_median = bg.bgmed;
  s.background_mean = [];
  s.bg_type = 'dark_on_light';
  s.n_bg_std_thresh_low = nan;
  s.istouching = td.istouching;
  
  varargout = cell(size(varargin));
  for i=1:nargout
    varargout{i} = s.(varargin{i});
  end
  return;
end
  
  
  
if nargin == 1,
  readall = true;
else
  readall = false;
  varargout = cell(1,nargin-1);
end;

fid = fopen(filename,'rb');
if fid < 0
  return
end

while true,

  s = fgetl(fid);
  if strcmp(s,'end header') || ~ischar(s),
    break;
  end

  [param,value] = read_line(s,fid);
  if isempty(param),
    continue;
  end;

  if readall,
    params.(param) = value;
  else
    varargout = set_output(param,value,varargin,varargout);
  end;

end;

if readall,
  varargout{1} = params;
end;

fclose(fid);



function out = set_output(param,value,in,out)

for i = 1:length(in),
  if strcmp(in{i},param),
    out{i} = value;
  end;
end;

function [param,value] = read_line(s,fid)

i = strfind(s,':');
if isempty(i),
  param = [];
  value = [];
  return;
end;
param = s(1:i-1);
value = s(i+1:end);

specialparams = {'background median','background mean',...
                 'background mad','background std',...
                 'hfnorm','fracframesisback',...
                 'background center','background dev',...
                 'isarena'};
stringparams = {'bg_algorithm','version','expbgfgmodel_filename','movie_name','data format','bg_type'};
pickledparams = {'roipolygons'};

isspecial = ismember(param,specialparams);
isstring = ismember(param,stringparams);
ispickled = ismember(param,pickledparams);

if isspecial,
  sz = str2double(value);
  value = fread(fid,sz/8,'double');
elseif ispickled,
  sz = str2double(value);
  value = fread(fid,sz,'char');  
elseif isstring,
  % leave value as string
else
  tmp = str2double(value);
  if ~isempty(tmp) && ~isnan( tmp ),
    value = tmp;
  end;
end;

param = strrep(param,' ','_');
