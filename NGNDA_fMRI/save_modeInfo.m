function save_modeInfo(filename,RS_modeInfo,modeinfo_loc)
% saves mode info file to correct location
save(filename,'RS_modeInfo');
movefile(filename,modeinfo_loc);

end