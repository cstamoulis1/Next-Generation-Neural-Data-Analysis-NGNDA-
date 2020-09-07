function master_submit_sc2(filenum, segstart,segend,ncores)
addpath('~/MIToolbox/')
addpath('~/ENSEMBLE_EMD')

%Submits one batch job using the cluster profile and the mpi queue
%sc2_masterjob parallelizes the segment-specific analysis

c=parcluster;
c.AdditionalProperties.WallTime = '48:00:00'; %to be changed by the user as needed
c.AdditionalProperties.QueueName = 'mpi'; %when more than 18 cores (CPUs), i.e., more than one node are requested, the mpi queue needs to be used 
  c.AdditionalProperties.AdditionalSubmitArgs = '--mem-per-cpu=2G' %the mode decomposition step is memory intensive, so at least 2G of memory should typically be requested; 
% more memory should be requested in the case of substantial spike rates in the data 

c.saveProfile

  l = filenum;
nc = ncores;
s1 = segstart;
s2 = segend;
  
c.batch(@sc2_masterjob,0,{l, s1,s2},'Pool',nc);

end
