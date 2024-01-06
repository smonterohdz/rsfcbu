function [newsubmask] = ImproveSeed_FC_2(xImg1,submask)
%ImproveSeed_FC_2 This is David's method for clustering
%   Detailed explanation goes here


%%
% CLUSTER STEP 1
% Get the image correlation matrix

C = xImg1;
R = corrcoef( C );


%%
% CLUSTER STEP 2
% Let's cluster my way
%
% NOTE THAT
% For each Iteration - corI(iIter): R2, Cluster and xImg1sub if I want to restart from that point
% I also need - corBase: xImg1
% for plotting I need - corBase: nVpial, lstV, fvB
% corBase.xImg1 = xImg1; corBase.fvB = fvB; corBase.lstV = lstV;
% corBase.nVpial = nVpial;
%
% save corBase corI
% ii=1;R2=corI(ii).R2; Cluster=corI(ii).Cluster; xImg1sub=corI(ii).xImg1sub;

R2 = R;

nV = size(xImg1,2);

Cluster = [1:nV]';

xImg1sub = xImg1;


%corBase.xImg1 = xImg1; 
%corBase.fvB = fvB; 
%corBase.lstV = lstV;
%corBase.nVpial = nVpial;

iIter = 0;
corI = [];


%%
% CLUSTER STEP 3
% GO FAST
% I speed things up to start with by sweeping through and
% pairing every vertix with its most correlated vertix. I am finding two
% sweeps of this reduces the number of clusters from >8000 to ~500, which
% is close to the number of measurements. At that point, we move on to the
% next stage. I even ran a third time.
for iIter=1:2
    R2(find(abs(R2-1)<1e-5)) = -2;
    for ir = 1:nV

        maxR2 = max(R2(ir,:));
        ic = find( abs(R2(ir,:)-maxR2)<1e-5, 1 );

        if Cluster(ir)~=Cluster(ic)
            iAll = find( Cluster==Cluster(ir) );
            iAll2 = find( Cluster==Cluster(ic) );

            iir2 = [iAll' iAll2'];
            fprintf('%d\n',ir) ;
            min_iir2 = min(iir2);
            Cluster(iir2) = min_iir2;

            xTmp = mean(xImg1(:,iir2),2);
            xImg1sub(:,iir2) = xTmp * ones(1,length(iir2));
        end
    end

    R2 = corrcoef( xImg1sub );

    %iIter = iIter + 1;
    corI(iIter).R2 = R2;
    corI(iIter).Cluster = Cluster;
    corI(iIter).xImg1sub = xImg1sub;
    corI(iIter).maxR2 = maxR2;
    corI(iIter).nUniqueClusters = length(unique(Cluster));

    %corI(iIter)
end

%%
% CLUSTER STEP 4
% Each step find cluster pair with max R2 

nSteps = 10;

for iLoop = 1:nSteps
    R2(find(abs(R2-1)<1e-5)) = -2;

    maxR2 = max(R2(:));
    
    [ir,ic] = find( abs(R2-maxR2)<1e-5, 1  );
    
    if Cluster(ir)==Cluster(ic)
        disp('keyboard');
        keyboard
        %finish the execution of the loop
    else
        iAll = find( Cluster==Cluster(ir) );
        iAll2 = find( Cluster==Cluster(ic) );
        
        iir2 = [iAll' iAll2'];
        fprintf('%d - %.5f; length(iir2)=%d\n',iLoop,maxR2,length(iir2) )
        min_iir2 = min(iir2);
        Cluster(iir2) = min_iir2;
        
        xTmp = mean(xImg1(:,iir2),2);
        xImg1sub(:,iir2) = xTmp * ones(1,length(iir2));
        
        R2tmp = mean( (xTmp-mean(xTmp))*ones(1,nV) .* (xImg1sub-ones(size(xTmp,1),1)*mean(xImg1sub,1)), 1) ./ ( std(xTmp,1) * std(xImg1sub,1,1) );
        
        for ii=1:length(iir2)
            R2(iir2(ii),:) = R2tmp;
            R2(:,iir2(ii)) = R2tmp';
        end
    end
    
end

iIter = iIter + 1;
corI(iIter).R2 = R2;
corI(iIter).Cluster = Cluster;
corI(iIter).xImg1sub = xImg1sub;
corI(iIter).maxR2 = maxR2;
corI(iIter).nUniqueClusters = length(unique(Cluster));

newsubmask = submask;
newsubmask.groups = Cluster;

%corI(iIter)


end