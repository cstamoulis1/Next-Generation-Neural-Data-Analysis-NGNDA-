function save_reduced_connectomes(filename,RS_net,connectome_loc_new)
% saves reduced connectome to current location and moves file to desired
% directory
save(filename,'RS_net');
movefile(filename,connectome_loc_new);
end