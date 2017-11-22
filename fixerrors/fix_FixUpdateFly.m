function fix_FixUpdateFly(handles,fly)
% sets fly plot properties based on fly data
% splintered from fixerrorsgui 6/21/12 JAB

if isdummytrk(handles.trx(fly))
  return;
end

fix_SetFlyVisible(handles,fly,'on');
ii = handles.trx(fly).off + handles.f;
if isalive(handles.trx(fly),handles.f)
  i = ii;
else
  i = 1;
end

x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
b = 2*handles.trx(fly).b(i);
theta = handles.trx(fly).theta(i);
ellipseupdate(handles.hellipse(fly),a,b,x,y,theta);
xwingl = handles.trx(fly).xwingl(i);
ywingl = handles.trx(fly).ywingl(i);
wingladj = 0.85*[xwingl-x ywingl-y]+[x y];
xwingladj = wingladj(1);
ywingladj = wingladj(2);
xwingr = handles.trx(fly).xwingr(i);
ywingr = handles.trx(fly).ywingr(i);
wingradj = 1.15*[xwingr-x ywingr-y]+[x y];
xwingradj = wingradj(1);
ywingradj = wingradj(2);

xleft = x - b*cos(theta+pi/2);
yleft = y - b*sin(theta+pi/2);
xright = x + b*cos(theta+pi/2);
yright = y + b*sin(theta+pi/2);
xhead = x + a*cos(theta);
yhead = y + a*sin(theta);
xtail = x - a*cos(theta);
ytail = y - a*sin(theta);

set(handles.htailmarker(fly),'xdata',[xtail,x],'ydata',[ytail,y]);
set(handles.hleft(fly),'xdata',xleft,'ydata',yleft);
set(handles.hright(fly),'xdata',xright,'ydata',yright);
set(handles.hhead(fly),'xdata',xhead,'ydata',yhead);
set(handles.htail(fly),'xdata',xtail,'ydata',ytail);
set(handles.hwingl(fly),'xdata',xwingladj,'ydata',ywingladj);
set(handles.hwingr(fly),'xdata',xwingradj,'ydata',ywingradj);
set(handles.hwinglinel(fly),'xdata',[x xwingladj],'ydata',[y ywingladj]);
set(handles.hwingliner(fly),'xdata',[x xwingradj],'ydata',[y ywingradj]);

i0 = ii - floor((handles.nframesplot-1)/2);
i1 = ii + handles.nframesplot - 1;
i0 = max(i0,1);
i1 = min(i1,handles.trx(fly).nframes);
set(handles.hpath(fly),'xdata',handles.trx(fly).x(i0:i1),...
  'ydata',handles.trx(fly).y(i0:i1));

% handles.needssaving = 1;
%guidata( handles.figure1, handles )

