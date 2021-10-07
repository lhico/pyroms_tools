addpath(genpath('/home/dsasaki/models/roms/20210901/matlab'))

% -- part A -- %
% the number of grid points must be divisible by ref_ration at each direction
ref_ratio=3;
Istr = 20;
Iend = 32;
Jstr = 98;
Jend = 111;
F=coarse2fine('bacia_santos_donor.nc',  'bacia_santos_rec1.nc', ref_ratio,Istr,Iend,Jstr,Jend);

Gnames={'bacia_santos_donor.nc',  'bacia_santos_rec1.nc'};
[S,G]=contact(Gnames,'roms_ngc.nc');

ref_ratio=5;
Istr = 10;
Iend = 25;
Jstr = 9;
Jend = 24;
F=coarse2fine('bacia_santos_rec1.nc',  'bacia_santos_rec2.nc', ref_ratio,Istr,Iend,Jstr,Jend);

$ -- part B -- %
% update grids as necessary and use them
Gnames={'bacia_santos_donord.nc',  'bacia_santos_ref1d.nc',  'bacia_santos_ref2c.nc'};
[S,G]=contact(Gnames,'roms_ngc_2.nc')