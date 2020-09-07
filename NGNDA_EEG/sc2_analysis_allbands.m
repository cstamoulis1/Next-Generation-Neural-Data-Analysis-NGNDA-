function sc2_analysis_allbands(file_num, seg_num)

%purpose: to analyze continuous neurophysiological data and compute connectivity using
%multiple methods (cross-correlation and Mutual Information (MI))
%Allows for true parallelization in time using high-resolution processing windows

%INPUTS: file_num - file number to be processed; seg_num; segment number to be processed

%OUTPUTS: saves the mode characteristics and mode-specific cross-correlation, lag, and MI matrices
 
% Authors: Yuna Lee and Caterina Stamoulis, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School
% Last version: December 31, 2018

%Calls the complete EEMD approach
% Torres, M. E., Colominas, M. A., Schlotthauer, G., and Flandrin, P. (2011, May):
% A complete ensemble empirical mode decomposition with adaptive noise
% (ICASSP), 2011 (pp. 4144-4147). DOI: 10.1109/ICASSP.2011.594726


load electrode_labels %variable name lab1; each patient may have different electrodes (particularly in invasive studies)

load sample_rate_PtX %variable sample_rate; vector of sampling rates for each file; sometimes the data acquisition sampling rate is switched mid-study
 
numfile = file_num;
fs = sample_rate(numfile);
npts = 2*fs; %2-s sliding analysis window
k = seg_num;


nchan = length(lab1); 

%Create pointers to all channels (softload the data and create slice_matrices
%Parallelization will be done for these matrices

%Create pointers to all channels
for i = 1:nchan
str1 = sprintf('filtdata_file%g_chan%s_new.mat',numfile, lab1{i});
pointfile = matfile(str1);
str2 = sprintf('pointer%g = pointfile;',i);
eval([str2])

%calculate globar statistics, including 99th percentile for artifacts
quant_global(i,:) = quantile(abs(pointfile.new_mat), [0.25 0.50 0.99]); 
end

nseg = floor(length(pointfile.new_mat)/npts);
q_upper = median(quant_global(:,3));				  
clear quant_global

for i = 1:nchan %loop to create matrix/amplitude thresholds
str3 = sprintf('sig1 = pointer%g.new_mat(1,(k-1)*npts+1:k*npts);',i);
eval([str3])
slice_mat(i,:) = sig1; %creates the k-th matrix slice
end


for i = 1:nchan
%Computer quantiles and extreme outliers to detect artifacts
%calculate local (signal-specific/segment-specific statistics)
q_sig = quantile(slice_mat(i,:), [0.25 0.50 0.99]);
if(q_sig(3) >q_upper)
  slice_mat(i,:) = NaN; %artifacts
end
end


nmodes = 10; %typically less than 10 modes have meaningful data
mode_mat=cell(nchan,nmodes,npts); %typically 
				  mode_freq=cell(nmodes,nchan);
				  mode_ampl=cell(nmodes,nchan);				  
				  mi_mode = cell(nmodes, nchan, nchan);
correl_mode = cell(nmodes, nchan, nchan);
lag_mode = cell(nmodes, nchan, nchan);

%Parallel processing at the channel level given the time it takes to do the mode decompositions

for i = 1:nchan
i
seg1 = slice_mat(i,:);
if(max(isnan(seg1))==1) %all vector entries are NaN
   md1 = nan(nmodes,npts); %typically < 10 modes have useful info
   else
     md1 = ceemdan(seg1,0.2, 200,1000);

 end
				  mode_mat{i}(:,:) = md1(1:nmodes,:);

   end
   

for m = 1:nmodes 
   for l = 1:nchan
				  seg1 = mode_mat{l}(m,:);
				  if(max(isnan(seg1))==0) 
				    mode_freq{m}(l) = (fs/2)*length(crossing(seg1))/npts;
				  mode_ampl{m}(l) = max(abs(seg1));
				  else
mode_freq{m}(l) = NaN;
				  mode_ampl{m}(l) = NaN;
				  end
				  for n = 1:nchan
				  seg2 = mode_mat{n}(m,:);
				  
				  if(max(isnan(seg1))==0 && max(isnan(seg2))==0)
				    mi_mode{m}(l,n) = mi(seg1',seg2');
[corr1 lag1] = xcorr(seg1,seg2,'coeff');
[mx1 pos1] = max(corr1);
				  correl_mode{m}(l,n) = mx1;
				  lag_mode{m}(l,n) = lag1(pos1);
 else
   mi_mode{m}(l,n) = NaN;
				  correl_mode{m}(l,n) = NaN;
				  lag_mode{m}(l,n) = NaN;
end %closes the if-else
end %closes the n-loop
end %closes the l-loop
end %closes the mode loop

res_struct1.freq = mode_freq;
res_struct1.ampl = mode_ampl;
res_struct1.MI = mi_mode;
res_struct1.correl = correl_mode;
res_struct1.lag = lag_mode;

clear mode_freq mode_ampl mi_mode correl_mode lab_mode slice_mat mode_mat
				  
str4 = sprintf('save struct_file%g_seg%g_sc2 res_struct1',numfile,k);
	       eval([str4])


				  end

