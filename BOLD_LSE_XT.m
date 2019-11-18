% Function to calculate difference of model and experimental BOLD data
% fit_parameters = {w_f, kappa, tau, alpha, beta}

% Created by: James Pang 
% Date of creation: Nov. 14, 2015

function difference = BOLD_LSE_XT(fit_parameters, astrocyte_model, ...
                                p, subject, calibration_index, extent)

%% Load experiment data: either raw or smooth
% Load data from singleSourceData folder
% Make sure that input for subject matches the case of data file name
% if strcmpi(data_type, 'raw')
%     load(cat(2,'singleSourceData/rawTS', subject, '.mat'))
%     ExperimentalResponse = avgTSER;
% elseif strcmpi(data_type, 'smooth')
%     load(cat(2,'singleSourceData/', subject, '.mat'))
%     ExperimentalResponse = smoothTSER;
% end

load(cat(2,'singleSourceData/', subject, '.mat'))
ExperimentalResponse = smoothTSER;



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
                        

%% Get experimental response at any x and t

ExperimentalResponse = ExperimentalResponse;
time_experiment = 0.25:0.25:20;     % 20 seconds with 250ms interval
space_experiment = distanceInterp*1e-3;  % need to multiply by 1e-3
[time_experiment_mat, space_experiment_mat] = meshgrid(time_experiment, space_experiment);

%% Interpolate predicted response to experimental time

time_predicted = PredictedResponse.time;
space_predicted = PredictedResponse.space;
[time_predicted_mat, space_predicted_mat] = meshgrid(time_predicted, space_predicted);

PredictedResponse = interp2(time_predicted_mat, space_predicted_mat, real(PredictedResponse.Y_xt), ...
    time_experiment_mat, space_experiment_mat);

%% Subtract mean to measure fluctuations around mean

PredictedResponse = PredictedResponse  - mean(PredictedResponse(:));


%% Normalize so that experimental and predicted response have the same height

ExperimentalResponse = ExperimentalResponse/max(ExperimentalResponse(:));
PredictedResponse = PredictedResponse/max(PredictedResponse(:));

%% First 5 values are not calibrated properly so take difference after n index
%%% Cut the measured difference from the spatial extent specified

[~, ind_left] = min(abs(space_experiment/1e-3 - extent(1))); 
[~, ind_right] = min(abs(space_experiment/1e-3 - extent(2))); 

n = calibration_index;
difference = PredictedResponse(ind_left:ind_right,n+1:end) - ExperimentalResponse(ind_left:ind_right,n+1:end);
difference = difference(:);
