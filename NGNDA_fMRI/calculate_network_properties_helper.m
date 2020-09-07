function props = calculate_network_properties_helper(small,cmat,thresholds,threshold_labels,measure)
    % purpose: To estimated network properties for each connectivity matrix for a set of desired thresholds
    % Calls multiple functions in the Sporns Brain Connectivity Toolbox
    % Brain Connectivity Toolbox: https://sites.google.com/site/bctnet/
    % Reference: Complex network measures of brain connectivity: Uses and interpretations.
    % Rubinov M, Sporns O (2010) NeuroImage 52:1059-69.
    % Author of this function: Sean Parks, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School 
    % Last version: September 5, 2020 
      % INPUT: small = 1 or 0 (yes or no) depending on whether user wants to estimate small worldness of the network (computationally expensive)
    %                thresholds input
    %        cmat = connectivity matrix 
    %        thresholds = threshold for connectivity matrix in order to westimate weighted adjacency matrices. In these matrices values below the threshold are set to zero
    %       Elements in the adjacency matrix at or above the threshold are the original connectivity values
    %        threshold_labels = cell array containing name of each
    %                           threshold, length of cell array is the same
    %                           as the length of the threshold array
    
    %        
    
    % OUTPUT: props = structure that contains network properties for
    %                 connectivity matrix for each threshold provided

    
    pathToSpornsCode = ''; % path to sporns code in order to use network properties functions
    % path is of form '/pathtoSpornsCode/SpornsCodeDirectory'
    addpath(pathToSpornsCode);
        
    cmat(isnan(cmat)) = 0; % any nans set to 0

    % Mutual information is a probabilistic measure. Therefore, mi(i,j) may be slightly different than mi(j,i), leading to a non-symmetric connectivity matrix.
      % Since some the network property estimations expect a symmetric connectivity and adjacency matrix, the upper triangular matrix of the origian MI matrix
      % is used to create a new symmetric MI matrix
    if measure == 3 % if MI connectivity matrix, must make symmetric due to stochasticity in MI calculation
        cmat_diag_vec = diag(cmat); % creates vector of main diagonal
        cmat_diag_mat = diag(cmat_diag_vec); % creates zero matrix except for main diagonal
        cmat_triu = triu(cmat,1); % takes upper triangular of matrix above diagonal
        cmat = cmat_triu + cmat_triu' + cmat_diag_mat; % reconstructs connectivity matrix by adding together upper triangular, its transpose, and the diagonal of the original connectivity matrix
    end
    
    
    A = ones(size(cmat));
    N = size(A,1);
    props = cell(1,length(thresholds)); % cell array that contains cell for each threshold input
    
    for thr = 1:length(thresholds) % goes through each threshold
        thresh = thresholds(thr); 
        A_t = cmat.*(cmat >= thresh); % thresholds connectivity matrix
        % any values above threshold remain, any below are set to 0
        A_t = weight_conversion(A_t,'normalize'); % normalizes the matrix by dividing by max connectivity value in matrix
        A_t = round(A_t,3); % rounds to correct for round off error so that matrix is symmetric
        props{thr} = struct; % creates empty structure in cell
        props{thr}.threshold = thresh; % stores threshold value in structure as well as threshold label in structure
        props{thr}.threshold_label = threshold_labels(thr);
        
        if any(A_t(logical(triu(A,1))) > 0) % if any values are nonzero (i.e., there are connections to analyze)

            props{thr}.mean_connectivity = nanmean(A_t(logical(triu(A,1))),'all'); % calculates mean of connection strengths in upper triangular (not including diagonal)
            sec1 = A_t(logical(triu(A,1)));
            props{thr}.median_connectivity = nanmedian(sec1(sec1 ~= 0),'all'); % computes median connection strength out of nonzero connectivities in upper triangular (not including diagonal)
        
            A_t_bin = weight_conversion(A_t,'binarize'); % creates binary version of adjacency matrix (all non zero set to 1)
            lambdas_t = svd(A_t_bin); 
            props{thr}.natural_connectivity = double(log(sum(exp(vpa(lambdas_t,5000)))/N)); % computes natural connectivity of matrix
	    % The natural connectivity quantifies the redundancy of alternative
            % routes in the network by evaluating the weighted number of closed walks of all lengths and can be seen as
            % an average eigenvalue obtained from the graph spectrum. 
	    % Jun Wu1, Mauricio Barahona, Yuejin Tan, Hongzhong Deng, Robustness of Regular Graphs Based on Natural Connectivity, https://www.researchgate.net/publication/46585739_Robustness_of_Random_Graphs_Based_on_Natural_Connectivity, 2010.
            
			
			
            rt_b = weight_conversion(A_t_bin,'autofix');  % removes self loops (puts 0s on diagonal) for binary matrix
        
            rand_clust = NaN(20,1); % compute 20 global clustering coefficients and characteristic path lengths from 20 random graphs from the binary connectivity matrix rt_b
            rand_path = NaN(20,1);
            
            if small % if small worldness needs be computed
                parfor iter = 1:20 % go through 20 iterations
        % random graph
                    addpath(pathToSpornsCode); % add sporns code path so each parallel worker can access the functions 
                    [rand_b, ~] = randmio_und(rt_b, 50); % generates random graph from binary connectivity matrix rt_b (randomizes the connections)
                    rand_b = weight_conversion(rand_b,'autofix'); % ensures 0s on the diagonal (necessary for Sporns network property codes)
                    clusterCoef_local_rand = clustering_coef_bu(rand_b); % computes local clustering coefficient of random graph
                    rand_clust(iter) = nanmean(clusterCoef_local_rand(~isinf(clusterCoef_local_rand)),'all'); % mean of local clustering coefficient across all nodes is global clustering coefficient
                    distances_rand = distance_bin(rand_b); % computes lengths of shortest paths between each pair of nodes
                    rand_path(iter) = charpath(distances_rand,0,0); % calculates characteristic path length of random graph
                end
                
                clusterCoef_global_rand = nanmean(rand_clust); % takes mean global clustering coefficient and mean path length across all random graphs generated
                characteristic_path_length_rand = nanmean(rand_path);

                clusterCoef_local_bin = clustering_coef_bu(rt_b); % computes global clustering coefficient of binary connectivity graph from earlier
                clusterCoef_global_bin = nanmean(clusterCoef_local_bin(~isinf(clusterCoef_local_bin)),'all'); 
                distances_R = distance_bin(rt_b);
                characteristic_path_length_bin = charpath(distances_R,0,0); % computes characteristic path length of binary connectivity graph from earlier 

        % small worldedness calculation
                gamma_R = clusterCoef_global_bin/clusterCoef_global_rand; % computes and stores small worldness of this connectivity matrix for this particular threshold
                lambda_R = characteristic_path_length_bin/characteristic_path_length_rand;
                props{thr}.small_worldness = gamma_R/lambda_R;
            end
            
            % computes transitivity, node degree, global efficiency, global
            % cluster coefficient, modularity, eigenvector centrality on
            % weighted undirected connectivity matrix with 0s on diagonal
            props{thr}.transitivity = transitivity_wu(weight_conversion(A_t,'autofix')); 
            props{thr}.degreeNodes = degrees_und(weight_conversion(A_t,'autofix'));
            props{thr}.globalEfficiency = efficiency_wei(weight_conversion(A_t,'autofix'));
            props{thr}.clusterCoef_local = clustering_coef_wu(weight_conversion(A_t,'autofix'));
            props{thr}.clusterCoef_global = nanmean(props{thr}.clusterCoef_local,'all');
            [~,Q_t]=modularity_und(weight_conversion(A_t,'autofix'));
            props{thr}.modularity = Q_t;
            props{thr}.centrality_eigenvector = eigenvector_centrality_und(weight_conversion(A_t,'autofix'));
            
            % computes first svd eigenvalue of weighted undirected
            % connectivity matrix with nonzero entries on diagonal
            svd_t = svd(A_t);
            props{thr}.svd_eigenvalue = svd_t(1);
        end   
    end
end
