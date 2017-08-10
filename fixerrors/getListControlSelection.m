function s = getListControlSelection(h)
assert(isa(h,'matlab.ui.control.UIControl'));
v = h.Value;
s = h.String{v};
