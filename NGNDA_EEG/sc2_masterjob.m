function sc2_masterjob(filenum, seg_start,seg_end)

%True parallelization in time (segment)
addpath('~/MIToolbox/')
addpath('~/ENSEMBLE_EMD')


  parfor k = seg_start:seg_end

  sc2_analysis_allbands(filenum, k);
end

  delete(gcp('nocreate'))

  end
  
