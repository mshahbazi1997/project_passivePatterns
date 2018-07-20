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
    case '5shades'
        colours                 = {'black','maroon','medred','orange','yellow'};
        canvas                  = 'blackonwhite';
        opt.general.markertype  = 'o';
        opt.general.markersize  = 6;
        opt.dot.facealpha       = 0.7;
        opt.dot.edgealpha       = 1;
        opt.dot.sizedata        = 40;
        opt.dot.gapwidth        = [1 1 1 1];
    case '5shadesNoMarkers'
        colours                 = {'black','maroon','medred','orange','yellow'};
        canvas                  = 'blackonwhite';
        opt.general.markertype  = 'none';  
    case '1black'
        %colours                 = {'black'};
        colours                 = plt.helper.get_shades(6,'gray','decrease',10);
        canvas                  = 'blackonwhite';
        %opt.general.markerfill = {[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0]};
        opt.general.markertype  = 'o';
        opt.general.markersize  = 6;
        opt.dot.facealpha       = 0.7;
        opt.dot.edgealpha       = 1;
        opt.dot.sizedata        = 40;
        opt.dot.gapwidth        = [1 0.75 1 1];

end;

