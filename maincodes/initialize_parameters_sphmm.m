function parameters = initialize_parameters_sphmm(data, parameters_of_initialize)
% data of shape (nS, nT, nC)
shape = size(data);
nC = shape(3);

num_subj = shape(1);
% get the settings of initialize step
num_cluster = parameters_of_initialize.num_cluster;
win_size = parameters_of_initialize.win_size;
if isempty(parameters_of_initialize.sample_size)
    sample_size = 0.5;
else
    sample_size = parameters_of_initialize.sample_size;
end

if isempty(parameters_of_initialize.method)
    method = 'sliding_window';
else
    method = parameters_of_initialize.method;  % "window"
end


if strcmp(method, 'sliding_window')
    disp('Applying sliding window intialization');
    % Init_cov = zeros(nC, nC, num_cluster);
    % Init_inv = Init_cov;
    Init_mu = zeros(nC, num_cluster);
    num_sample = ceil(sample_size*num_subj);
    fprintf('%.0f subjects are chosen for initialize\n', num_sample);
    % win_size = 30;
    %lambda = 0.2;
    % use sliding window to get the cov and Init_mu
    
    temp_perm = randperm(num_subj);
    sample_idx = temp_perm(1:num_sample);
    sample_data = data(sample_idx,:,:);
    sliding_method = 'L1';  % None
    if parameters_of_initialize.allcov
        [OutputCov, Index, all_cov] = sliding_window(sample_data, num_cluster, win_size, sliding_method);
    else
        [OutputCov, Index] = sliding_window(sample_data, num_cluster, win_size, sliding_method);
    end
    
    parameters.initcov = permute(OutputCov, [2,3,1]);
    numWin = size(Index, 2);
    TCmean = zeros(num_sample, numWin, nC);
    for isubj = 1: num_sample
        window_range = 1:30;
        for j=1: numWin
            TCmean(isubj, j,:) = mean(sample_data(isubj, window_range,:),2);
            window_range = window_range + 1;
        end
    end
    
    mean_vec = reshape(TCmean, [num_sample*numWin, nC]);
    Index = reshape(Index, [num_sample*numWin,1]);
    if parameters_of_initialize.allcov
        parameters.index = Index;
        parameters.allcov = all_cov;
    end
    for icluster = 1: num_cluster
        i_index = (Index == icluster);
        if sum(i_index)==1  % only one timecourse was marked within one cluster
            Init_mu(:, icluster) = mean_vec(i_index,:);
        else
            Init_mu(:, icluster) = mean(mean_vec(i_index,:), 1);
        end
    end
    parameters.initMu = Init_mu;
elseif strcmp(method, 'window')
    
else
    error('Please specify the initialize method...');
    
    
    
end
% calculate the mean according to the Index from slding window


%[Emp_corr,slide_mu] = slide_init(sample_data,num_cluster,win_size,overlap);
% initialize the cov and empirical in lasso
%     glasso_cov = Emp_corr;
%     glasso_theta = Emp_corr;
%     for j=1:num_cluster
%         lambda = sqrt(2*log(num_roi))/(num_time/num_cluster);
%         %lambda = 0.1;
%         [glasso_cov(:,:,j), glasso_theta(:,:,j)] = graphlasso(Emp_corr(:,:,j), lambda);
%     end
%     Init_mu{i} = slide_mu';
%     Init_cov{i} = glasso_cov;
%     Init_inv{i} = glasso_theta;
% end


end