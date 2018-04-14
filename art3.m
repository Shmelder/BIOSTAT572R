function ART_out = art3(Y,X, ARTplot, constant, gamma, B, dB)
    % X: a n*p matrix, where p is the number of variables and n is the
    %     sample size (no intercept);
    % Y: a n*1 vector;
    % ARTplot = 1 if display the histogram for ART; and 0, o.w. The default
    %     value is 0.
    % consant: a vector of choices of constants considered in lambda.
    %          The default value is constant = 0.5:0.05:4; 
    % gamma: the significance level. default value gamma = 0.05;
    % B: the number of bootstrap samples. Default value B=1000;
    % dB: the number of double bootstrap samples. Default value dB=1000;
    % ART_out: a two element row vector, where the first element is the 
    %     index of the selected covariate, and the second element is te 
    %     ART p-value.
 
    
    if nargout > 2
        error('myfuns:art:TooManyOutputs', ...
            'requires at most 2 onputs');
    end
    
    if nargin > 7
        error('myfuns:art:TooManyInputs', ...
            'requires at most 2 optional inputs');
    end
    
    
    % setting default values
    default = {0 0.5:0.05:4 0.05 1000 1000};
    if nargin == 2
        [ARTplot, constant, gamma, B, dB] = default{:};
    elseif nargin == 3
        [constant, gamma, B, dB] = default{2:5};
    elseif nargin == 4
        [gamma, B, dB] = default{3:5};
    elseif nargin == 5
        [B, dB] = default{4:5};
    elseif nargin == 6
        dB = default{5};
    end
    if isempty(ARTplot)==1
        ARTplot = default{1};
    end
    if isempty(constant)==1
        CPBplot = default{2};
    end
    if isempty(gamma)==1
        gamma = default{3};
    end
    if isempty(B)==1
        B = default{4};
    end
    
    %% AE : This code chunk runs a linear regression for every covariate
    %% and finds the maximum of these covariates (most significant predictor)

    [n,p] = size(X);
    beta = zeros(2,p); err = zeros(1,p);
    for k = 1:p
        beta(:,k) = regress(Y, [ones(n,1) X(:,k)]);
        err(k) = mean((Y - [ones(n,1) X(:,k)]*beta(:,k)).^2);
    end
    
    %% AE : This code finds the t-statistic and checks it against
    %% the lambda value as part of the procedure that estimates
    %% the limiting distribution

    [temp,I] = sort(err);
    khat = I(1);
    betahat = beta(:,khat);
    ehat = Y - [ones(n,1) X(:,khat)]*betahat;
    v1hat = (sum(ehat.^2)/(n-2))/(mean(X(:,khat).^2)-(mean(X(:,khat)))^2);
    tnT = abs(sqrt(n)*betahat(2)/sqrt(v1hat));
    
    lambda = max(constant*sqrt(log(n)), norminv(1-gamma/(2*p),0,1)); 
    cpb = zeros(B,1); 
    GridIndexT = zeros(B,length(constant)); 
    Vgrid = zeros(B, 1);
        
reject = zeros(1,length(constant));
for bs = 1: B
    bs
    bsidx = randsample(n, n, true);
    bsY = Y(bsidx);
    bsX = X(bsidx,:);
            
    bsRSS = var(bsY,1)-(mean(bsY)*mean(bsX,1)- ...
        mean(bsxfun(@times, bsX, bsY),1)).^2./var(bsX,1,1);
    [bsRSShat, bskhat] = min(bsRSS);
    bsbetahat = regress(bsY, [ones(n,1) bsX(:,bskhat)]);
            
    %%%%%%%% CPB
    cpb(bs) = sqrt(n)*(bsbetahat(2) - betahat(2));
    
    %%%% ART
    bsehat = bsY - [ones(n,1) bsX(:,bskhat)]*bsbetahat;
            
    bsv1hat = (sum(bsehat.^2)/(n-2))/...
        (mean(bsX(:,bskhat).^2)-(mean(bsX(:,bskhat)))^2);
    bstnT = abs(sqrt(n)*bsbetahat(2)/sqrt(bsv1hat));
        
    GridIndexT(bs,:) = (max(bstnT,tnT) <= lambda);
            
    %bsCov = (bsX'*bsX)/n - (mean(bsX))'*mean(bsX);
    bsCov=cov(bsX,1);
    bsVar = diag(bsCov);
                
    bsZ = (sqrt(n)*(mean(bsxfun(@times, ehat(bsidx), bsX),1) ...
        - mean(ehat(bsidx))*mean(bsX,1)...
        - (mean(bsxfun(@times,ehat,X),1)-mean(ehat)*mean(bsX,1))))';
                
    temp = limitFun(bsZ, bsVar, bsCov, 0, p);
    Vgrid(bs)=temp(1);
        
    dcpb = zeros(dB,1);
    dGridIndexT = zeros(dB,length(constant)); 
    dVgrid = zeros(dB, 1);
    
    %% AE : This code chunk uses a double bootstrap to calculate the 
    %% tuning parameter lambda used in the limiting distribution.      

    for dbs=1:dB
        dbsidx = randsample(n, n, true);
        dbsY = bsY(dbsidx);
        dbsX = bsX(dbsidx,:);
            
        dbsRSS = var(dbsY,1)-(mean(dbsY)*mean(dbsX,1)- ...
            mean(bsxfun(@times, dbsX, dbsY),1)).^2./var(dbsX,1,1);
        [dbsRSShat, dbskhat] = min(dbsRSS);
        dbsbetahat = regress(dbsY, [ones(n,1) dbsX(:,dbskhat)]);
            
            %%%%%%%% CPB
        dcpb(dbs) = sqrt(n)*(dbsbetahat(2) - bsbetahat(2));
        
            %%%% ART
        dbsehat = dbsY - [ones(n,1) dbsX(:,dbskhat)]*dbsbetahat;
            
        dbsv1hat = (sum(dbsehat.^2)/(n-2))/...
            (mean(dbsX(:,dbskhat).^2)-(mean(dbsX(:,dbskhat)))^2);
        dbstnT = abs(sqrt(n)*dbsbetahat(2)/sqrt(dbsv1hat));
        
        dGridIndexT(dbs,:) = (max(dbstnT,bstnT) <= lambda);
               
             %bsCov = (bsX'*bsX)/n - (mean(bsX))'*mean(bsX);
        dbsCov=cov(dbsX,1);
        dbsVar = diag(dbsCov);
                 
        dbsZ = (sqrt(n)*(mean(bsxfun(@times, bsehat(dbsidx), dbsX),1) ...
            - mean(bsehat(dbsidx))*mean(dbsX,1)...
            - (mean(bsxfun(@times,bsehat,bsX),1)-mean(bsehat)*mean(dbsX,1))))';
        
        temp = limitFun(dbsZ, dbsVar, dbsCov, 0, p);
        dVgrid(dbs)=temp(1);
    end
       
    for j=1:length(constant)
        dCIrange = bsxfun(@plus, dcpb.*(1-dGridIndexT(:,j)), ...
            bsxfun(@times, dVgrid, dGridIndexT(:,j)));
        sorttemp = sort(dCIrange,1);
        lowerquantile = min(sorttemp(dB*gamma/2,:));
        upperquantile = max(sorttemp(dB*(1-gamma/2),:));
    
        if ((sqrt(n)*(bsbetahat(2)-betahat(2))>upperquantile)||...
            (sqrt(n)*(bsbetahat(2)-betahat(2))<lowerquantile))
            reject(j) = reject(j)+1;
        end
               
    end
end
reject = reject/dB;
[temp,optJ] = min(abs(bsxfun(@minus, reject, gamma)),[],2);
      
       
CIrange = bsxfun(@plus, cpb.*(1-GridIndexT(:,optJ)), ...
    bsxfun(@times, Vgrid, GridIndexT(:,optJ)));
ART_pvalue = 2*min(sum(CIrange <= sqrt(n)*betahat(2)), ...
           sum(CIrange >= sqrt(n)*betahat(2)))/B;
    
    if nargout == 0
        fprintf('The ART p-value is %0.4f.\n',ART_pvalue)
    else
        ART_out = [khat ART_pvalue];
    end
    
    if (ARTplot ==1)
        hist(CIrange,20)
        temp = hist(CIrange,20);
        hold on
        plot(sqrt(n)*betahat(2)*ones(2,1), [0 max(temp)+1],'-.')
        hold off
        title('ART','fontsize',13)
        ylabel('Frequency','fontsize',12)
        xlabel('$A_n^*$','fontsize',12, 'Interpreter','latex')
    end
end

function bsV = limitFun(Z, variance, covariance, b, p)

    V = bsxfun(@plus, Z./variance, b*(bsxfun(@rdivide, ...
        covariance, variance)-1));
    M = bsxfun(@rdivide, bsxfun(@plus, Z, b*covariance).^2, variance);

    [temp, I] = sort(M);
    maxK = I(p,:);
    bsV = zeros(1,p);
    for j = 1: p
        bsV(j) = V(maxK(j),j);
    end
end


