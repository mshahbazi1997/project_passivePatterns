function [canvas,colours,opt] = pp1_style(styID)
%% Description
%   Default plotting parameters specified by the style styID
%   Define additional styles within the switch statement. Users can create
%   their own styles for their projects and use the style package to point
%   towards the user-defined style
%
%
% Author
%   Naveed Ejaz (ejaz.naveed@gmail.com)

canvas           = 'blackonwhite';
opt              = [];
opt.save.journal = 'brain';

switch(styID)
    case 'default'
        colours                 = {'black','lightgray'};
        opt.display.ax          = 'normal';
    case '5fingers'
        colours                 = {[0.2 0.14 0.53] [0.29 0.49 0.36] [0.97 0 0] [0 0.58 0.77] [0.8 0.27 0.5]};
        opt.display.ax          = 'normal';
        opt.hist.facealpha      = 0.4;
        opt.general.markersize  = 6;
        opt.general.linewidth   = 1.3;
    case 'numDigits'
        %colours                 = {[0 0 0] [0.5 0 0] [0.9 0 0] [1 0.6 0] [0 0 0]};
        colours                 = {[0 0 0] [0.5 0 0] [0.9 0 0] [1 0.6 0] [0.4 0.4 0.4]};
        opt.display.ax          = 'normal';
        opt.hist.facealpha      = 0.4;
        opt.general.markersize  = 6;
        opt.general.linewidth   = 1.3;
end;

