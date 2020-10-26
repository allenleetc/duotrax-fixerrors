function s = getListControlSelection(h)

% Added by Allen Lee for DTFE in 2017

assert(isa(h,'matlab.ui.control.UIControl'));
v = h.Value;
s = h.String{v};
