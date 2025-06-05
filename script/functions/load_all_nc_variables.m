function [data_struct,info_struct] = load_all_nc_variables(ncfile)
% [data_struct,info_struct] = load_all_nc_variables(ncfile);
%%
% ncfile = '/Users/natbrig/Dropbox/School/Matlab/BiogeochemicalArgo/data/bodc/3901496/3901496_meta.nc';

info_struct = ncinfo(ncfile);

for vix = 1:length(info_struct.Variables)
    vname = info_struct.Variables(vix).Name;
    vname2 = vname;
    vname2(vname==' ') = '_';
    data_struct.(vname2) = ncread(ncfile,vname);
end
