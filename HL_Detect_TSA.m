function HL_est = HL_Detect_TSA(vid_name,N,nBins,Zscore)

    obj = VideoReader(vid_name);
    vid = read(obj);
    
    %% Get required parameters
    RGB = vid(:,:,:,1);
    % Get number of frames in the input video as nof
    nof=obj.NumFrames;
    % Get height (h) and width (w) of the frame
    [h,w,~]=size(RGB);

    % Set number of maximum iterations (M) to 3.
    M = 3;
    
    x_ref=w/2;
    [Xmat_or,Ymat_or]=meshgrid(1:w,1:h);
    t = 1:w;
    BW_or = zeros(h,w);
    delta_h = 0.025*h;

    minh = round(h*0.075);

    %% Listener Block
    % Set frame index to 1
    k=1;
    % Call HLDA to get state [y_k, theta_k] of the HL at k=1. Append the
    % current state to HL estimates matrix HL_est.
    HL_est(k,:)=HLDA(RGB);
    % Calculate the error metric
    err(k,:) = getError(HL_est(k,:),RGB,Xmat_or,Ymat_or,nBins,x_ref,delta_h);
    % After the 1st frame use ROI (referenced in XXX) and repeat the process for k = 2 to N 
    for k = 2:N
        % Get the kth frame
        RGB = vid(:,:,:,k);
        % Determine the ROI
        [ROI,h_lims] = ROI_RECT(RGB,HL_est(k-1,:),minh,x_ref,h);
        % Call HLDA to get state
        x = HLDA(ROI);
        % Update the vertical position w.r.t. image coordinate frame
        x(1)=x(1)+h_lims(1);
        % Append the current state to HL estimates matrix HL_est.
        HL_est(k,:)=x;
        % Calculate the error metric
        err(k,:) = getError(HL_est(k,:),RGB,Xmat_or,Ymat_or,nBins,x_ref,delta_h);
        % Plot the estimated HL.
        HLplot(RGB,HL_est(k,:),t,x_ref);
    end
    
    %% Create univariate autoregressive integrated moving average (ARIMA) models.
    % For y, ARIMA(2,2,0) and GARCH(1,1)
    Mdly = arima('ARLags',1:2,'D',2,'Variance',garch(1,1));
    % For theta, ARIMA(2,0,0) and GARCH(1,1)
    Mdlt = arima('ARLags',1:2,'D',0,'Variance',garch(1,1));

    % Optimize model coefficients
    y_Mdl_est = estimate(Mdly,HL_est(end-N+Mdly.P+1:end,1),'Display','off');
    theta_Mdl_est = estimate(Mdlt,HL_est(end-N+Mdlt.P+1:end,2),'Display','off');

    % Set AD (Absence Detector variable) false
    AD = 0;

    %% Main Block   
    while k<nof
        k=k+1;
        % Get the frame
        RGB = vid(:,:,:,k);
        % Check absence of HL
        if AD>0
            delh = round(minh);
            ROI = RGB(h-2*delh:h,:,:);
            BW = BW_or; BW(h-2*delh:h,:)=1;
    
            x = HLDA(ROI);
            HL_est(k,:) = [x(1)+h-2*delh x(2)];
            err(k,:) = getError(HL_est(k,:),RGB,Xmat_or,Ymat_or,nBins,x_ref,delta_h);
            % abs((err(k,:)-median(err(k-10:k-1,:)))./median(err(k-10:k-1,:)))
            abs((err(k,:)-err(indTemp))./err(indTemp))
            if abs((err(k,:)-err(indTemp))./err(indTemp))>0.05
                HL_est(k,:)=[h,0];
                HLplot(RGB,HL_est(k,:),t,x_ref);
                ROIplot(BW);
                drawnow;
                continue;
            else
                AD = 0;
                HLplot(RGB,HL_est(k,:),t,x_ref);
                ROIplot(BW);
                drawnow;
                continue;
            end
        end
        
        % TS Block
        % Get forecasts for y and theta (yf, tf) and their variances
        % (yf_var, tf_var)
        [yf,yf_var] = forecast(y_Mdl_est,1,HL_est(end-N+Mdly.P+1:end,1));
        [tf,tf_var] = forecast(theta_Mdl_est,1,HL_est(end-N+Mdlt.P+1:end,2));
        % Calculate the standard deviations
        yf_sd = sqrt(yf_var); tf_sd = sqrt(tf_var);

        % Determine the ROI
        [ROI,AD,BW,delh] = ROI_TSM(RGB,[yf,tf],yf_sd,tf_sd,Zscore,w);
        if AD<1
            % Call HLDA to get state
            x = HLDA(ROI);
            % Append the updated state w.r.t. image coordinate system to HL
            % estimates matrix HL_est.
            HL_est(k,:) = [yf-delh+x(1)/cosd(tf), tf+x(2)];
            % Calculate the error metric
            err(k,:) = getError(HL_est(k,:),RGB,Xmat_or,Ymat_or,nBins,x_ref,delta_h);
            
            
            % Control Loop
            % Set no. of iterations (NOI) to 0.
            NOI = 0;
            minh2 = 0; 
            
            while NOI < M && (abs((err(k,:)-mean(err(end-4:end,:)))./mean(err(end-4:end,:)))>0.15 || x(1)>3*h || x(1)<0)
                NOI = NOI + 1;

                [ROI,h_lims] = ROI_RECT(RGB,HL_est(k-1,:),minh2,x_ref,h);%HL_est(i,:)
                BW = BW_or; BW(h_lims(1):h_lims(2),:)=1;
                delh = (-h_lims(1)+h_lims(2))/2;
               
                x = HLDA(ROI);
                HL_est(k,:) = [x(1)+h_lims(1) x(2)];

                err(k,:) = getError(HL_est(k,:),RGB,Xmat_or,Ymat_or,nBins,x_ref,delta_h);
                minh2=minh2+round(0.075*h);
            end
            HLplot(RGB,HL_est(k,:),t,x_ref);
            ROIplot(BW);
            drawnow;
            
            
        else
            AD = 1;
            HL_est(k,:) = [h, 0]; %% aşağı yukarı karar verecek bişiler yaz.
            indTemp = k-1;
            HLplot(RGB,HL_est(k,:),t,x_ref);
            BW = BW_or; BW(h-2*round(minh):h,:)=1;
            ROIplot(BW);
            drawnow;
        end
    end
end

function ERR=getError(x,RGB,Xmat_or,Ymat_or,nBins,x_ref,delta_h)
    [ind_sky,ind_sea]=find_indices_of_regions(x,x_ref,Xmat_or,Ymat_or,delta_h);
    Im = rgb2gray(RGB);
    ERR = -sqrt(mean((histcounts(Im(ind_sky),nBins,'Normalization','probability')-histcounts(Im(ind_sea),nBins,'Normalization','probability')).^2));
    if isnan(ERR)
        ERR = 0;
    end
end

function [ind_sky,ind_sea]=find_indices_of_regions(x,x_ref,Xmat,Ymat,delta_h)
    yI = x(1); theta =x(2);
    c = yI-tand(theta)*x_ref;
    ind_ref = tand(theta)*Xmat+c;

    ind_lower = tand(theta)*Xmat+c-delta_h;
    ind_upper = tand(theta)*Xmat+c+delta_h;

    ind_sky = Ymat<ind_ref & Ymat>ind_lower;
    ind_sea = Ymat>ind_ref & Ymat<ind_upper;
end


%% Horizon Line Detection Algorithm
function x = HLDA(RGB)
    [~,threshOut] = edge(rgb2gray(RGB),"canny");
    thrNew = threshOut*2;
    thrNew(thrNew>=1)=0.9999;
    BW = edge(rgb2gray(RGB),"canny",thrNew);
    BW2 = bwareaopen(BW,round(sum(BW(:))*0.005));

    theta = 0:0.25:180-0.25;
    [Rad,xp] = radon(BW2,theta);
    
    [row_peak,col_peak] = find(ismember(Rad,max(max(Rad))));
    dist = xp(row_peak);
    th = theta(col_peak);

    dist = dist(1); th = th(1);
    
    x(2) = 90-th;
        if dist~=0
            uy = sign(dist)*sind(th);
            A = abs(dist)/uy;
        else
            A = 0;
        end
    x(1) = -A+size(RGB,1)/2;
end

%% ROI algorithms
function [ROI,AD,BW,delh] = ROI_TSM(RGB,x,yf_sd,tf_sd,Zscore,w)
    AD = 0;
    HL = tand(x(2))*[1 w]+x(1)-tand(x(2))*w/2;
    HL_angled = tand(x(2)+tf_sd)*w+x(1)-tand(x(2)+tf_sd)*w/2;
    HL = round(HL);

    delh = round(Zscore*(yf_sd + abs(HL_angled-HL(2))));
    r = [HL(1)-delh HL(1)+delh HL(2)+delh HL(2)-delh];

    c = [0 0 w w];
    BW = roipoly(RGB,c,r);
    S = sum(BW);

    if sum(S)/prod(size(BW))>0.0025

        % width cropping
        ind = S<max(S(:));
        BW(:,ind)=0;
    
        R = RGB(:,:,1); 
        G = RGB(:,:,2);
        B = RGB(:,:,3);
    
        ROI(:,:,1) = reshape(R(BW),[],sum(~ind));
        ROI(:,:,2) = reshape(G(BW),[],sum(~ind));
        ROI(:,:,3) = reshape(B(BW),[],sum(~ind));
    else
        ROI = []; AD = 1; BW = []; delh = 0;
    end
end

function [ROI,h_lims] = ROI_RECT(RGB,x,minh,x_ref,h)
    y=round(tand(x(2))*(1:2*x_ref)+x(1)-tand(x(2))*x_ref);
    h_lims = [min(y)-minh max(y)+minh];
    h_lims(h_lims>h)=h; h_lims(h_lims<1)=1;
    ROI = RGB(h_lims(1):h_lims(2),:,:);
end
%% plotting functions
function ROIplot(BW)
    B = bwboundaries(BW);
    bound = B{1,1};
    hold on;
    plot(bound(:,2),bound(:,1),'y',LineWidth=3);
    hold off;
end

function HLplot(RGB,x,t,x_ref)
    imshow(RGB); 
    hold on;
    y=round(tand(x(2))*t+x(1)-tand(x(2))*x_ref);
    line([1,t(end)],[y(1),y(end)],'Color','r','LineWidth',3);
    drawnow;
    hold off;
end