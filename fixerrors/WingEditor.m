classdef WingEditor < handle
  % Drag to adjust wing angles on two flies in an axis
  
  properties
    hFig
    hAx
    hPts % [2x2] graphics handles. hPts(fly,iwing) (iwing=1/2 for l/r)
    hLines % [2x2] etc. hLines(fly,iWing)
    
    xyFly % [2x2] xyFly(fly,icoord) (icoord=1/2 for x/y)
    aFly % [2x1] aFly(fly)
    thetaFly % [2x1] thetaFly(fly)
        
    tfActivated
    dragFly % 0 (no drag in progress), 1, or 2
    dragWing % 1 or 2
    wingtipConAFac = [4*.85 4*1.15]; % wingtipConAFac(fly)
  end
  
  methods
    function obj = WingEditor(hpts,hlines)
      szassert(hpts,[2 2]);
      szassert(hlines,[2 2]);
      
      obj.hFig = ancestor(hpts(1),'figure');
      obj.hAx = ancestor(hpts(1),'axes');
      
      obj.hPts = hpts;
      obj.hLines = hlines;
      obj.xyFly = nan(2,2);
      obj.aFly = nan(2,1);
      obj.thetaFly = nan(2,1);
            
      obj.tfActivated = false;
      obj.dragFly = 0;
      obj.dragWing = nan;      
    end
  end
  
  methods
        
    function angAbs = getWingAbsAngles(obj)
      % angAbs: [2x2] absolute angles (relative to fixed coord, not 
      % relative to fly). andAbs(fly,iwing)

      hpts = obj.hPts;
      xyfly = obj.xyFly;
      angAbs = nan(2,2);
      for fly=1:2
      for iwing=1:2
        h = hpts(fly,iwing);
        angAbs(fly,iwing) = atan2(get(h,'YData')-xyfly(fly,2),...
                                  get(h,'XData')-xyfly(fly,1));
      end
      end
    end
    
    function collapse(obj,fly)
      xy = obj.xyFly(fly,:);
      a = obj.aFly(fly,:);
      th = obj.thetaFly(fly,:);
      thTail = th+pi;
      xyTail = xy + 4*a*[cos(thTail) sin(thTail)];

      for iwing=1:2
        afac = obj.wingtipConAFac(iwing);
        xyWingtip = wingtipConFcn(xyTail,xy,a,afac);
        set(obj.hPts(fly,iwing),'XData',xyWingtip(1),'YData',xyWingtip(2));
        
        hline = obj.hLines(fly,iwing);
        hlinexdata = get(hline,'XData');
        hlineydata = get(hline,'YData');
        hlinexdata(2) = xyWingtip(1);
        hlineydata(2) = xyWingtip(2);
        set(hline,'XData',hlinexdata,'YData',hlineydata);
      end
    end
    
    function ptBDF(obj,src,evt)
      assert(obj.tfActivated);
      [fly,iwing] = find(src==obj.hPts);
      if isempty(fly)
        [fly,iwing] = find(src==obj.hLines);
      end
      obj.dragFly = fly;
      obj.dragWing = iwing;
    end
    
    function figWBMF(obj,src,evt)
      fly = obj.dragFly;
      iwing = obj.dragWing;
      if fly>0
        xyCur = obj.hAx.CurrentPoint(1,:);
        xyCur = xyCur(1:2);
        afac = obj.wingtipConAFac(iwing);
        xy = wingtipConFcn(xyCur,obj.xyFly(fly,:),obj.aFly(fly),afac);
        set(obj.hPts(fly,iwing),'XData',xy(1),'YData',xy(2));
        
        hline = obj.hLines(fly,iwing);
        hlinexdata = get(hline,'XData');
        hlineydata = get(hline,'YData');
        hlinexdata(2) = xy(1);
        hlineydata(2) = xy(2);
        set(hline,'XData',hlinexdata,'YData',hlineydata);
      end
    end
    
    function figWBUF(obj,src,evt)
      obj.dragFly = 0;
      obj.dragWing = nan;
    end
    
  end
  
end
