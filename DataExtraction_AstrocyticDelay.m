%%% Code to extract all relevant model and experimental data for each subject

function data = DataExtraction_AstrocyticDelay(subject, calibration_index)

%% Initialization and loading of mat files
p = params;     %% Load model parameters

if calibration_index==0
    filename = ['DATA+FIGURES_K_wf=0.1to1_DisregardFirstFiveDataPoints',...
        '/DATA_',subject,'_AstrocyticDelay_IncludeFirstFivePoints.mat'];
elseif calibration_index==5
    filename = ['DATA+FIGURES_K_wf=0.1to1_DisregardFirstFiveDataPoints',...
        '/DATA_',subject,'_AstrocyticDelay_DisregardFirstFivePoints.mat'];
end

load(filename);     %% Loading prediction/model results

load singleSourceData/editedFileParams.mat      %% Loading v_b and Gamma parameters

load(cat(2,'singleSourceData/rawTS', subject, '.mat'))  %% Loading experiment results
load(cat(2,'singleSourceData/', subject, '.mat'))

%% Extracting v_b and Gamma parameter values
data.v_b = Vmat(find(strcmp(associatedParams, subject)));      % Note that 1e-3 is needed for future usage
data.Gamma = -Gmat(find(strcmp(associatedParams, subject))); 

%% Extracting raw and smooth experimental data
data.experiment.time = 0.25:0.25:20;     % 20 seconds with 250ms interval
data.experiment.space = distanceInterp;
data.experiment.BOLDraw_all = avgTSER;
data.experiment.BOLDsmooth_all = smoothTSER;
data.experiment.BOLDraw_center = avgTSER(8,:);
data.experiment.BOLDsmooth_center = smoothTSER(8,:);

%% Extracting predicted response and calculations
data.prediction.models = models;
for n_models = 1:8
    if strcmp(data.prediction.models{n_models}.data_type, 'raw')
        maxHeight_experiment_center = max(data.experiment.BOLDraw_center);  %% Taking maximum of experiment BOLD
    elseif strcmp(data.prediction.models{n_models}.data_type, 'smooth')
        maxHeight_experiment_center = max(data.experiment.BOLDsmooth_center);
    end
    
    prediction_parameters = [data.prediction.models{n_models}.parameters.wf, ...
                             data.prediction.models{n_models}.parameters.kappa, ...
                             data.prediction.models{n_models}.parameters.tau];
    PredictedResponse = HemodynamicResponse1D(prediction_parameters, ...
                        data.prediction.models{n_models}.astrocytic_model, p, ...
                        data.v_b*1e-3, data.Gamma);       %% Model result
    
    data.prediction.time = PredictedResponse.time;
    data.prediction.space = PredictedResponse.space;

    centerResponse_prediction = interp1(data.prediction.time, ...
                real(PredictedResponse.Y_xt(size(PredictedResponse.Y_xt,1)/2+1,:)), ...
                data.experiment.time);      % interpolate predicted response to experimental time
    centerResponse_prediction = centerResponse_prediction - mean(centerResponse_prediction);    % subtract to measure fluctuations around mean       
    centerResponse_prediction = (centerResponse_prediction/max(centerResponse_prediction))*maxHeight_experiment_center;   %normalize so that experimental and predicted response have the same height
    
    data.prediction.models{n_models}.BOLD_all = real(PredictedResponse.Y_xt); %y: distance, x:time
    data.prediction.models{n_models}.BOLD_center = centerResponse_prediction;
end


