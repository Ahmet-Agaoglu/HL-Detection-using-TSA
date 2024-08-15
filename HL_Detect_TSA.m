function HL_est = HL_Detect_TSA(vid_name,N)
    % HL_Detect_TSA: Detects the horizon line in a video using time series analysis.
    % Inputs:
    %   vid_name - Name of the video file.
    %   N - Number of frames to process in the Listener Block.
    % Outputs:
    %   HL_est - Matrix of estimated horizon line states. Rows correspond
    %   to the frame indeces. The 1st colomn is the vertical position of the
    %   Horizon Line (HL) in pixels whereas the 2nd column consists the orientation
    %   angle in angles.

    % Read the video file.
    obj = VideoReader(vid_name);
    vid = read(obj);
    
    %% Extract video parameters
    RGB = vid(:,:,:,1); % Extract the first frame.
    nof = obj.NumFrames; % Get the number of frames in the video.
    [h,w,~] = size(RGB); % Get the height and width of the frame.

    % Set algorithm parameters
    M = 3; % Maximum number of iterations for the control loop.
    nBins = 0:255; % Histogram bins vector for error calculation.
    Zscore = 1.96; % Standard z-score corresponding to 95% confidence interval.
    
    % Define coordinate matrices for the pixels.
    [Xmat_or,Ymat_or] = meshgrid(1:w,1:h);
    
    % Initialize a blank binary image for ROI calculations.
    BW_or = zeros(h,w);
    
    % Set the height for the parallelogram used in error calculation.
    delta_h = 0.025 * h;
    
    % Set the minimum height for the rectangular ROI.
    minh = round(h * 0.075);

    %% Listener Block
    [HL_est, err, k] = ListenerBlock(vid, N, Xmat_or, Ymat_or, nBins, w, delta_h, minh, h);

    %% Create univariate autoregressive integrated moving average (ARIMA) models.
    % For y, ARIMA(2,2,0) and GARCH(1,1)
    Mdly = arima('ARLags',1:2,'D',2,'Variance',garch(1,1));
    % For theta, ARIMA(2,0,0) and GARCH(1,1)
    Mdlt = arima('ARLags',1:2,'D',0,'Variance',garch(1,1));

    %% Optimize model coefficients
    y_Mdl_est = estimate(Mdly,HL_est(end-N+Mdly.P+1:end,1),'Display','off');
    theta_Mdl_est = estimate(Mdlt,HL_est(end-N+Mdlt.P+1:end,2),'Display','off');

    % Set AD (Absence Detector variable) false
    AD = 0;
    
    %% Process the remaining frames
    while k<nof
        k=k+1;
        % Get the frame
        RGB = vid(:,:,:,k);
        
        if AD>0 % Check absence of HL
            % if HL is absent, route to Presence Detector Block
            [HL_est, err, AD] = PresenceDetector(RGB, HL_est, err, minh, BW_or, k, indTemp, Xmat_or, Ymat_or, nBins, w, h, delta_h);
            continue;
        end
        % if HL is present, route to Main Block
        %% Main Block
        % Initialize the TS Block
        % Get forecasts for y and theta (yf, tf) and their variances (yf_var, tf_var)
        [yf, yf_var] = forecast(y_Mdl_est, 1, HL_est(end-N+Mdly.P+1:end, 1));
        [tf, tf_var] = forecast(theta_Mdl_est, 1, HL_est(end-N+Mdlt.P+1:end, 2));
        
        % Calculate the standard deviations
        yf_sd = sqrt(yf_var); 
        tf_sd = sqrt(tf_var);
        
        % Determine the ROI based on the forecasted states and their variances
        [ROI, AD, BW, delh] = ROI_TSM(RGB, [yf, tf], yf_sd, tf_sd, Zscore, w);
        
        if AD < 1 % Control Loop
            % Call HLDA to get the estimated state of the HL within the ROI
            x = HLDA(ROI);
            
            % Update the state estimate relative to the image coordinate system and 
            % append it to the HL estimates matrix HL_est
            HL_est(k, :) = [yf - delh + x(1)/cosd(tf), tf + x(2)];
            
            % Calculate the error metric based on the current state estimate
            err(k, :) = getError(HL_est(k, :), RGB, Xmat_or, Ymat_or, nBins, w, delta_h);
            
            % Initialize Control Loop variables
            NOI = 0; % Number of iterations
            minh2 = 0; % Rectangular ROI height increment
            
            % Control Loop: Iterate while error is significant or state estimates are out of bounds
            while NOI < M && (abs((err(k,:) - mean(err(end-4:end,:))) ./ mean(err(end-4:end,:))) > 0.15 || x(1) > 3*h || x(1) < 0)
                NOI = NOI + 1;
        
                % Adjust ROI to a larger rectangular region
                [ROI, h_lims] = ROI_RECT(RGB, HL_est(k-1, :), minh2, w, h);
                
                % Update the binary image BW for the new ROI
                BW = BW_or; 
                BW(h_lims(1):h_lims(2), :) = 1;
                % delh = (-h_lims(1) + h_lims(2)) / 2;
                
                % Re-estimate the HL state within the updated ROI
                x = HLDA(ROI);
                HL_est(k, :) = [x(1) + h_lims(1), x(2)];
                
                % Recalculate the error metric
                err(k, :) = getError(HL_est(k, :), RGB, Xmat_or, Ymat_or, nBins, w, delta_h);
                
                % Increment the ROI height adjustment for the next iteration
                minh2 = minh2 + round(0.075 * h);
            end
            
            % Plot the estimated HL and the ROI
            HLplot(RGB, HL_est(k, :), w);
            ROIplot(BW);
            drawnow;
            
        else % Absence Detector
            % Handle cases where the horizon line is not detected
            indTemp = k-1; % Record the last frame index when the HL was present.
            % Set the horizon line estimate based on the last known state
            if HL_est(k-1, 1) > h/2
                HL_est(k, :) = [h, 0];
                BW = BW_or; 
                BW(h-2*minh:h, :) = 1;
            else
                HL_est(k, :) = [1, 0];
                BW = BW_or; 
                BW(1:2*minh, :) = 1;
            end
            
            % Plot the fallback HL and the ROI
            HLplot(RGB, HL_est(k, :), w);            
            ROIplot(BW);
            drawnow;
        end
    end
end


    function [HL_est, err,k] = ListenerBlock(vid, N, Xmat_or, Ymat_or, nBins, w, delta_h, minh, h)
        % ListenerBlock: Processes the first N frames to initialize the system.
        % Inputs:
        %   vid - Video data.
        %   N - Number of frames to process.
        %   Other inputs are video parameters and initial settings.
        % Outputs:
        %   HL_est - Matrix of estimated horizon line states.
        %   err - Error metric for each frame.

        % Initialize frame index.
        k = 1;
        
        % Process the first frame using the full frame.
        RGB = vid(:,:,:,k); % Get the first frame.
        HL_est(k,:) = HLDA(RGB); % Estimate the horizon line state.
        err(k,:) = getError(HL_est(k,:), RGB, Xmat_or, Ymat_or, nBins, w, delta_h); % Calculate error.
        
        % Process the subsequent frames using a rectangular ROI.
        for k = 2:N
            RGB = vid(:,:,:,k); % Get the kth frame.
            [ROI, h_lims] = ROI_RECT(RGB, HL_est(k-1,:), minh, w, h); % Determine the ROI.
            x = HLDA(ROI); % Estimate the horizon line state within the ROI.
            x(1) = x(1) + h_lims(1); % Adjust the vertical position.
            HL_est(k,:) = x; % Append the current state to the estimates.
            err(k,:) = getError(HL_est(k,:), RGB, Xmat_or, Ymat_or, nBins, w, delta_h); % Calculate error.
            HLplot(RGB, HL_est(k,:), w); % Plot the estimated horizon line.
        end
    end

    function [HL_est, err, AD] = PresenceDetector(RGB, HL_est, err, minh, BW_or, k, indTemp, Xmat_or, Ymat_or, nBins, w, h, delta_h)
        % PresenceDetector: Detects the presence of the horizon line and adjusts the state estimates.
        % Inputs:
        %   RGB - The current video frame in RGB format.
        %   HL_est - The matrix containing the estimated states of the horizon line.
        %   err - The error metrics associated with the horizon line estimates.
        %   minh - Minimum height used for ROI calculation.
        %   BW_or - A blank binary image for ROI calculations.
        %   k - The current frame index.
        %   indTemp - The index of the last frame when the HL was present.
        %   Xmat_or, Ymat_or - Coordinate matrices of the pixels.
        %   nBins - Histogram bins vector used in error calculation.
        %   w, h - The width, height of the video frame.
        %   delta_h - The height of the parallelogram used in error calculation.
        % Outputs:
        %   HL_est - Updated matrix containing the estimated states of the horizon line.
        %   err - Updated error metrics associated with the horizon line estimates.
        %   AD - Absence Detector flag, indicating whether the horizon line was detected (0) or not (1).
    
        % Adjust the ROI based on the previous horizon line estimate
        if HL_est(k-1, 1) > 0
            % If the previous estimate is in the lower half of the frame, focus on the bottom region
            ROI = RGB(h - 2*minh:h, :, :);
            BW = BW_or; BW(h - 2*minh:h, :) = 1;
            x = HLDA(ROI);
            HL_est(k, :) = [x(1) + h - 2*minh, x(2)];
        else
            % Otherwise, focus on the top region
            ROI = RGB(1:2*minh, :, :);
            BW = BW_or; BW(1:2*minh, :) = 1;
            x = HLDA(ROI);
            HL_est(k, :) = x;
        end
    
        % Calculate the error metric for the current horizon line estimate
        err(k, :) = getError(HL_est(k, :), RGB, Xmat_or, Ymat_or, nBins, w, delta_h);
    
        % Check if the error exceeds a threshold compared to the reference frame
        if abs((err(k, :) - err(indTemp)) / err(indTemp)) > 0.05
            % If the error is too high, revert to the previous state estimate
            AD = 1;
            HL_est(k, :) = HL_est(k-1, :);
            HLplot(RGB, HL_est(k, :), w);
            ROIplot(BW);
            drawnow;
        else
            % Otherwise, update the Absence Detector flag and plot the current state
            AD = 0;
            HLplot(RGB, HL_est(k, :), w);
            ROIplot(BW);
            drawnow;
        end
    end

    function HLplot(RGB, x, w)
        % HLplot: Plots the estimated horizon line on the current video frame.
        % Inputs:
        %   RGB - The current video frame in RGB format.
        %   x - The estimated state of the horizon line [y_k, theta_k].
        %   w - The width of the video frame.
        
        % Display the current frame.
        imshow(RGB); 
        hold on;
        
        % Calculate the y-coordinates of the horizon line across the frame width.
        y = round(tand(x(2)) * (1:w) + x(1) - tand(x(2)) * w / 2);
        
        % Plot the horizon line in red with a thickness of 3 pixels.
        line([1, w], [y(1), y(end)], 'Color', 'r', 'LineWidth', 3);
        
        % Refresh the display to show the updated frame with the horizon line.
        drawnow;
        
        % Release the hold on the current figure.
        hold off;
    end

    function x = HLDA(RGB)
        % HLDA: Detects the horizon line in a given video frame or a ROI.
        % Inputs:
        % RGB - The current video frame in RGB format.
        % Outputs:
        %   x - The estimated state of the horizon line [y_k, theta_k]
        
        % Convert the RGB image to grayscale and perform Canny edge detection
        [~, threshOut] = edge(rgb2gray(RGB), "canny");
        
        % Set the threshold to eliminate the weak edges
        thrNew = threshOut * 2;
        thrNew(thrNew >= 1) = 0.9999;
        
        % Apply Canny edge detection with the new threshold
        BW = edge(rgb2gray(RGB), "canny", thrNew);
        
        % Remove small objects from the binary image
        BW2 = bwareaopen(BW, round(sum(BW(:)) * 0.005));
        
        % Define the angles for Radon transform
        theta = 0:0.25:180 - 0.25;
        
        % Perform the Radon transform on the binary image
        [Rad, xp] = radon(BW2, theta);
        
        % Find the peak in the Radon transform which corresponds to the HL
        [row_peak, col_peak] = find(ismember(Rad, max(max(Rad))));
        
        % Get the distance and angle corresponding to the peak
        dist = xp(row_peak);
        th = theta(col_peak);
        
        % If multiple peaks, choose the first one
        dist = dist(1); 
        th = th(1);
        
        % Calculate the orientation of the horizon line (theta_k)
        x(2) = 90 - th;
        
        % Calculate the vertical position of the horizon line (y_k)
        if dist ~= 0
            uy = sign(dist) * sind(th);
            A = abs(dist) / uy;
        else
            A = 0;
        end
        x(1) = -A + size(RGB, 1) / 2;
    end

    function ERR = getError(x, RGB, Xmat_or, Ymat_or, nBins, w, delta_h)
        % getError: Calculates the error metric between the sky and sea regions
        %           in the image based on the grayscale histograms.
        % Inputs:
        %   x - The estimated state of the horizon line [y_k, theta_k].
        %   RGB - The current video frame in RGB format.
        %   Xmat_or - X-coordinate matrix for the frame.
        %   Ymat_or - Y-coordinate matrix for the frame.
        %   nBins - Vector specifying histogram bins.
        %   w - The width of the video frame.
        %   delta_h - Height of the parallelogram used in error calculation.
        % Outputs:
        %   ERR - The error metric calculated between the sky and sea regions.
        
        % Identify the sky and sea regions based on the estimated HL state
        [ind_sky, ind_sea] = find_indices_of_regions(x, w, Xmat_or, Ymat_or, delta_h);
        
        % Convert the RGB image to grayscale
        Im = rgb2gray(RGB);
        
        % Calculate the histograms of the sky and sea regions
        hist_sky = histcounts(Im(ind_sky), nBins, 'Normalization', 'probability');
        hist_sea = histcounts(Im(ind_sea), nBins, 'Normalization', 'probability');
        
        % Compute the error as the negative square root of the mean squared difference
        ERR = -sqrt(mean((hist_sky - hist_sea).^2));
        
        % Handle cases where there is no sky or sea region.
        if isnan(ERR)
            ERR = 0;
        end
    end

    function [ind_sky, ind_sea] = find_indices_of_regions(x, w, Xmat, Ymat, delta_h)
        % find_indices_of_regions: Identifies the pixel indices for the sky and sea regions based on the estimated horizon line.
        % Inputs:
        %   x - The estimated state of the horizon line [y_k, theta_k].
        %   w - The width of the video frame.
        %   Xmat - X-coordinate matrix for the frame.
        %   Ymat - Y-coordinate matrix for the frame.
        %   delta_h - Height of the parallelogram used to define the region around the horizon line.
        % Outputs:
        %   ind_sky - Logical matrix indicating pixels belonging to the sky region.
        %   ind_sea - Logical matrix indicating pixels belonging to the sea region.
    
        % Calculate the reference line (horizon line) based on the estimated state.
        yI = x(1); 
        theta = x(2);
        c = yI - tand(theta) * w / 2;
        ind_ref = tand(theta) * Xmat + c;
    
        % Calculate the upper and lower bounds around the horizon line.
        ind_lower = tand(theta) * Xmat + c - delta_h;
        ind_upper = tand(theta) * Xmat + c + delta_h;
    
        % Identify the sky region as pixels above the reference line but below the upper bound.
        ind_sky = Ymat < ind_ref & Ymat > ind_lower;
        
        % Identify the sea region as pixels below the reference line but above the lower bound.
        ind_sea = Ymat > ind_ref & Ymat < ind_upper;
    end

    function [ROI, h_lims] = ROI_RECT(RGB, x, minh, w, h)
        % ROI_RECT: Determines a rectangular Region of Interest (ROI) around the estimated horizon line.
        % Inputs:
        %   RGB - The current video frame in RGB format.
        %   x - The estimated state of the horizon line [y_k, theta_k].
        %   minh - The minimum height added above and below the estimated horizon line.
        %   w - The width of the video frame.
        %   h - The height of the video frame.
        % Outputs:
        %   ROI - The cropped image corresponding to the Region of Interest.
        %   h_lims - The height limits defining the upper and lower boundaries of the ROI.
        
        % Calculate the vertical positions along the horizon line across the frame width.
        y = round(tand(x(2)) * (1:2 * w / 2) + x(1) - tand(x(2)) * w / 2);
        
        % Define the height limits of the ROI by extending minh pixels above and below the horizon line.
        h_lims = [min(y) - minh, max(y) + minh];
        
        % Ensure the height limits do not exceed the frame boundaries.
        h_lims(h_lims > h) = h; 
        h_lims(h_lims < 1) = 1;
        
        % Crop the image to create the ROI based on the height limits.
        ROI = RGB(h_lims(1):h_lims(2), :, :);
    end

    function [ROI, AD, BW, delh] = ROI_TSM(RGB, x, yf_sd, tf_sd, Zscore, w)
        % ROI_TSM: Generates the parallelogram ROI using Time Series Model forecasts.
        % Inputs:
        %   RGB - The current video frame in RGB format.
        %   x - The estimated state of the horizon line [y_k, theta_k].
        %   yf_sd - Standard deviation of the forecasted y position.
        %   tf_sd - Standard deviation of the forecasted theta (orientation).
        %   Zscore - Z-score value for the confidence interval (e.g., 1.96 for 95%).
        %   w - The width of the video frame.
        % Outputs:
        %   ROI - The extracted region of interest.
        %   AD - Absence Detector flag (1 if ROI is outside the frame, 0 otherwise).
        %   BW - Binary mask representing the ROI.
        %   Half-height of the ROI, dynamically calculated based on the input parameters.
    
        % Initialize Absence Detector (AD) flag to 0 (indicating the horizon is found).
        AD = 0;
        
        % Calculate the horizon line (HL) using the forecasted theta and y values.
        HL = tand(x(2)) * [1 w] + x(1) - tand(x(2)) * w / 2;
        
        % Calculate the horizon line for an angled forecast (theta + tf_sd).
        HL_angled = tand(x(2) + tf_sd) * w + x(1) - tand(x(2) + tf_sd) * w / 2;
        
        % Round the horizon line positions to integer values.
        HL = round(HL);
    
        % Calculate the vertical distance (delh) for the ROI height.
        delh = round(Zscore * (yf_sd + abs(HL_angled - HL(2))));
        
        % Define the ROI boundaries as a parallelogram.
        r = [HL(1) - delh, HL(1) + delh, HL(2) + delh, HL(2) - delh];
        c = [0 0 w w];
        BW = roipoly(RGB, c, r);
        
        % Calculate the sum of active pixels in the binary mask (BW).
        S = sum(BW);
    
        % Check if the ROI is within the frame.
        if sum(S) / numel(BW) > 0.0025
            % Crop the ROI based on the width of the active region.
            ind = S < max(S(:));
            BW(:, ind) = 0;
    
            % Extract the ROI from the RGB image using the binary mask.
            R = RGB(:, :, 1);
            G = RGB(:, :, 2);
            B = RGB(:, :, 3);
    
            ROI(:, :, 1) = reshape(R(BW), [], sum(~ind));
            ROI(:, :, 2) = reshape(G(BW), [], sum(~ind));
            ROI(:, :, 3) = reshape(B(BW), [], sum(~ind));
        else
            % If the ROI is outside the frame, set the Absence Detector flag to 1.
            ROI = [];
            AD = 1;
            BW = [];
            delh = 0;
        end
    end

    function ROIplot(BW)
        % ROIplot: Plots the Region of Interest (ROI) on the current video frame.
        % Inputs:
        %   BW - A binary image where active cells determine the ROI boundary.

        % Identify the boundaries of the ROI in the binary image
        B = bwboundaries(BW);

        % Extract the boundary coordinates
        bound = B{1,1};

        % Plot the ROI boundary on the current video frame
        hold on;
        plot(bound(:,2),bound(:,1),'y',LineWidth=3);
        hold off;
    end

