function xynew = wingtipConFcn(xy,xyFly,aFly,afac)
  % xynew must land on a circle of radius 4*aFly centered at xyFly
    
th = atan2(xy(2)-xyFly(2),xy(1)-xyFly(1));
xynew = afac*aFly*[cos(th) sin(th)]+xyFly;
