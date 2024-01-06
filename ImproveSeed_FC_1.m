function [newsubmask] = ImproveSeed_FC_1(seedHbX_ts,submask,r_thresh)
%ImproveSeed_1 This functions applies clustering to our submasks data
% version one utilizes Matlab clustering

ccoefseedHb = corrcoef(seedHbX_ts);
%remove the main diagonal with ones
ccoefseedHb = ccoefseedHb - eye(size(ccoefseedHb));

% and convert to a vector (as pdist)
dissimilarity = 1 - squareform(ccoefseedHb);

% decide on a cutoff
% remember that 0.4 corresponds to corr of 0.6!
%r_thresh = 0.8;

groups = 1;
while (length(unique(groups))<2)
    cutoff = r_thresh*2;

    %# perform complete linkage clustering
    Z = linkage(dissimilarity,'complete');

    % group the data into clusters
    % (cutoff is at a correlation of 0.5)
    groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
    if (length(unique(groups))<2) %if there's only one group, deccrease r_threh in 10 percent
        r_thresh = r_thresh*0.9;
    end
end

newsubmask = submask;
newsubmask.groups = groups;


% max_mu = 0;
% max_i = -1;
% 
% for i=1:length(unique(groups))
%     %obtain the submatrix with vertices within cluster i
%     ccoefseedHb_ = ccoefseedHb(groups==i,groups==i);
%     n_vert = size(ccoefseedHb_,1);
%     %obtain the lower triangular matrix without the diagonal
%     ccoefseedHb_ = ccoefseedHb_(find(tril(ones(n_vert),-1)));
%     % average the correlation values of vertices within group i
%     mu = mean(ccoefseedHb_);
%     if mu > max_mu
%         max_mu = mu;
%         max_i = i;
%         sum_submask_corr = sum(ccoefseedHb_);
%     end
% end
% idx_cc = groups==max_i;
% submask.mask_subsetseed = idx_cc;
% seed_HbX_ts = mean(seedHbX(:,idx_cc),2);

end