% h = myscatter(x,y,scale,color,[ncolors],[cm],['filled'],marker,...)
%
%  MYSCATTER Scatter/bubble plot.
%     MYSCATTER(X,Y,S,C) displays colored circles at the locations specified
%     by the vectors X and Y (which must be the same size).  This function
%     is similar to the SCATTER function. Instead of plotting objects as
%     patches, it instead plots them as lines, with a different line for
%     each color instead of a different patch for each point. The colors of
%     the lines will not change when the colormap of the axes is changed.
%     Instead, the colormap must be specified by the input. Some of the
%     differences between MYSCATTER and SCATTER are within the **'s below.
%  
%     S determines the area of each marker (in points^2). S can be a
%     vector the same length a X and Y or a scalar. If S is a scalar, 
%     MATLAB draws all the markers the same size. If S is empty, the
%     default size *12* is used.
%     
%     C determines the colors of the markers. When C is a vector the
%     same length as X and Y, the values in C are linearly mapped
%     to the colors in the *jet* colormap. When C is a 
%     length(X)-by-3 matrix, it directly specifies the colors of the  
%     markers as RGB values. C can also be a color string. See ColorSpec.
%  
%     MYSCATTER(X,Y) draws the markers in the default size and color.
%     MYSCATTER(X,Y,S) draws the markers at the specified sizes (S)
%     with a single color. This type of graph is also known as
%     a bubble plot.
%     *MYSCATTER(X,Y,S,C,NCOLORS) uses NCOLORS different colors in the
%     scatter.*
%     *MYSCATTER(X,Y,S,C,NCOLORS,COLORMAP) specifies the colors to use with
%     the COLORMAP input. If COLORMAP is a function handle, then NCOLORS
%     colors are generated by calling COLORMAP(NCOLORS). Otherwise,
%     COLORMAP must be a NCOLORS x 3 matrix, where COLORMAP(i,:) is the i
%     th color to use.*
%     MYSCATTER(...,M) uses the marker M instead of 'o'.
%     MYSCATTER(...,'filled') fills the markers.
%     *All extra arguments are fed into the plot command, and should come in
%     pairs. *
%  
%     H = MYSCATTER(...) returns handles to the scatter objects created.
%  
%     Use PLOT for single color, single marker size scatter plots.
%  
%     Example
%       load seamount
%       scatter(x,y,5,z)
function updatemyscatter(h,centers,x,y,s,c)

ncenters = length(centers);
[~,~,idx] = myhist(c,centers);
for i = 1:ncenters,
  set(h(i),'xdata',x(idx==i),'ydata',y(idx==i));
end

