function handles = fix_StorePanelPositions(handles)
% store the positions of the panels in the GUI
% splintered from fixerrorsgui 6/23/12 JAB

figpos = get(handles.figure1,'Position');

handles.upperrightpanel_tags = {'editpanel' 'swappanel'};
ntags = numel(handles.upperrightpanel_tags);
handles.upperrightpanel_dright = nan(1,ntags);
hnadles.upperrightpanel_dtop = nan(1,ntags);
for fni = 1:ntags
  fn = handles.upperrightpanel_tags{fni};
  h = handles.(fn);
  pos = get(h,'Position');
  handles.upperrightpanel_dright(fni) = figpos(3)-pos(1);
  handles.upperrightpanel_dtop(fni) = figpos(4)-pos(2);
end

handles.lowerrightpanel_tags = {'navigationpanel'};
ntags = numel(handles.lowerrightpanel_tags);
handles.lowerrightpanel_dright = nan(1,ntags);
for fni = 1:ntags
  fn = handles.lowerrightpanel_tags{fni};
  h = handles.(fn);
  pos = get(h,'Position');
  handles.lowerrightpanel_dright(fni) = figpos(3)-pos(1);
end
  
handles.bottom_tags = {'frameinfopanel' 'pnlNav'};
% ntags = numel(handles.bottom_tags);
% handles.bottom_width_norm = nan(1,ntags);
% handles.bottom_dleft_norm = nan(1,ntags);
% for fni = 1:ntags,
%   fn = handles.bottom_tags{fni};
%   h = handles.(fn);
%   pos = get(h,'Position');
%   handles.bottom_width_norm(fni) = pos(3)/figpos(3);
%   handles.bottom_dleft_norm(fni) = pos(1)/figpos(3);
% end

% axes: width expands to fit up against rightpanel_dright
% height expands, starting from existing y-pos

% pos = get(handles.mainaxes,'Position');
% % sliderpos = get(handles.frameslider,'Position');
% rightpanelpos = get(handles.seqinfopanel,'Position');
% handles.axes_dtop = figpos(4) - (pos(2)+pos(4));
% % handles.axes_dslider = pos(2) - (sliderpos(2)+sliderpos(4));
% handles.axes_drightpanels = rightpanelpos(1)-(pos(1)+pos(3));
