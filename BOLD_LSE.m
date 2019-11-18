% Function to calculate difference of model and experimental BOLD data
% fit_parameters = {w_f, kappa, tau, alpha, beta}

% Created by: James Pang 
% Date of creation: Nov. 14, 2015

function difference = BOLD_LSE(fit_parameters, astrocyte_model, ...
                                p, subject, data_type, calibration_index, plotting)

%% Load experiment data: either raw or smooth
% Load data from singleSourceData folder
% Make sure that input for subject matches the case of data file name
if strcmpi(data_type, 'raw')
    load(cat(2,'singleSourceData/rawTS', subject, '.mat'))
    ExperimentalResponse.Yt = avgTSER;
elseif strcmpi(data_type, 'smooth')
    load(cat(2,'singleSourceData/', subject, '.mat'))
    ExperimentalResponse.Yt = smoothTSER;
end



%% Load fileParams.mat that contains the estimated v_b and Gamma for subject

% load singleSourceData/fileParams.mat
load singleSourceData/editedFileParams.mat

% Find index of subject in associatedParams cell
index = find(strcmp(associatedParams, subject));

% Assigns the estimated experiment v_b and Gamma to variables for future use
v_b = Vmat(index);
Gamma = -Gmat(index);



%% Calculate the predicted response using the fitting parameters and 
%  experimental v_b and Gamma

PredictedResponse = HemodynamicResponse1D(fit_parameters, ...
                            astrocyte_model, p, v_b*1e-3, Gamma);
                        

%% Get experimental response at x==0

centerResponse_experiment = ExperimentalResponse.Yt(8,:);
maxHeight_experiment = max(centerResponse_experiment);
time_experiment = 0.25:0.25:20;     % 20 seconds with 250ms interval


%% Interpolate predicted response to experimental time

centerResponse_prediction = interp1(PredictedResponse.time, ...
    real(PredictedResponse.Y_xt(size(PredictedResponse.Y_xt,1)/2+1,:)), ...
    time_experiment);


%% Subtract mean to measure fluctuations around mean

centerResponse_prediction = centerResponse_prediction - mean(centerResponse_prediction);


%% Normalize so that experimental and predicted response have the same height

centerResponse_prediction = (centerResponse_prediction/max(centerResponse_prediction))*maxHeight_experiment;


%% First 5 values are not calibrated properly so take difference after n index

n = calibration_index;
difference = centerResponse_prediction(n+1:end) - centerResponse_experiment(n+1:end);


%% For plotting
if plotting
    RMSE = sqrt(sum(difference.^2));
    fig = figure('Visible','off');
%     subplot(1,2,1)
    plot(time_experiment(1:end), centerResponse_prediction(1:end), 'ks-', 'MarkerSize',5, 'MarkerFaceColor', 'k');
    hold on
    plot(time_experiment(1:end), centerResponse_experiment(1:end), 'ro-', 'MarkerSize',5, 'MarkerFaceColor', 'r');
    xlabel('time (s)','fontSize',18);ylabel('BOLD response (arb units)','fontSize',18);
    set(gca,'fontSize',18);
%     ylim([-1,1]);
%     title(['RMSE = ', num2str(RMSE)])
%     title('Comparison between predicted and experiment');
    legend('prediction', 'experiment')
    hold off
    
    set(fig, 'PaperPositionMode','auto')     % WYSIWYG
    %%%% MUST CHANGE FOR ALPHA AND BETA
    if calibration_index==0
        print('-depsc', [data_type,'_',subject,'_',astrocyte_model,'_wf='...
            num2str(fit_parameters(1)),'_kappa=',num2str(fit_parameters(2)),...
            '_tau=',num2str(fit_parameters(3)),'_RMSE=',num2str(RMSE),...
            '_IncludeFirstFivePoints.eps'])
    elseif calibration_index==5
        print('-depsc', [data_type,'_',subject,'_',astrocyte_model,'_wf='...
            num2str(fit_parameters(1)),'_kappa=',num2str(fit_parameters(2)),...
            '_tau=',num2str(fit_parameters(3)),'_RMSE=',num2str(RMSE),...
            '_DisregardFirstFivePoints.eps'])
    end
    
%     subplot(1,2,2)
%     plot(time_experiment(6:end), difference, 'b.-');
%     xlabel('time (s)','fontSize',18);ylabel('BOLD (arb units)','fontSize',18);
%     set(gca,'fontSize',18);title('Comparisons between predicted and experiment');
end
