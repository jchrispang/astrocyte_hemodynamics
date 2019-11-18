% Function to calculate 1D hemodynamic response for astrocytic simulations
% Created by: James Pang 
% Date of creation: Nov. 12, 2015
% Edited: May 30, 2016

function PredictedResponse = HemodynamicResponse1D(astrocyte_parameters, ...
                            astrocyte_model, p, v_b, Gamma)


%% NOTES
% Make sure that p = params
% Model parameters are in p
% Astrocyte_parameters is a vector of at most five parameters
% first: w_f
% second: kappa
% third: tau_d --> leave as [] if using astrocyte_model = 'original' 
% fourth: dt --> leave as [] if using astrocyte_model = 'original' and 'delay'

% Storing parameters on output data
PredictedResponse.v_b = v_b; 
PredictedResponse.Gamma = Gamma;


%% Setting up Fourier transform parameters and generate x and t vectors

% predefine vectors for temporal frequency
jvec = 0:p.Nw-1;
f = (jvec - p.Nw/2)/p.Nw*p.freqMax*2;
w = 2*pi*f;                                 % w = \omega

jveck = 0:p.Nk-1;
fk = (jveck - p.Nk/2)/p.Nk*p.spatialFreqMax*2;
k = 2*pi*fk;                                 % kx

t = 1/(f(2)-f(1))*1/p.Nw*(jvec - p.Nw/2);           % time vector
x = 1/(fk(2)-fk(1))*1/p.Nk*(jveck - p.Nk/2);        % position vector

[tt, xx] = meshgrid(t, x);          % mesh matrix of space and time
[wmat, kmat] = meshgrid(w, k);      % mesh matrix of spatial and temporal freq


%% Astrocytic model 

% Storing type of astrocytic model on output data
PredictedResponse.astrocyte_model = astrocyte_model;

switch astrocyte_model
    case 'no-delay'         % original model
        PredictedResponse.w_f = astrocyte_parameters(1);
        PredictedResponse.kappa = astrocyte_parameters(2);
        PredictedResponse.tau_d = [];  % storing parameters on output data
       
        H_t = 1;                        % astrocyte convolving function         
        
    case 'delay'            % astrocytic delay      
        PredictedResponse.w_f = astrocyte_parameters(1);
        PredictedResponse.kappa = astrocyte_parameters(2);
        PredictedResponse.tau_d = astrocyte_parameters(3);  % storing parameters on output data
        
        % Find index of t == tau_d 
        % Generalized for cases when t == tau_d is not exactly present
        tau_d = PredictedResponse.tau_d;
        [~, ind] = min(abs(t - tau_d));
        
        H_t = 1*(t==t(ind));                        % astrocyte convolving function 
    case 'low-pass_filter1'  % low-pass filter te^-t/dt      
        PredictedResponse.w_f = astrocyte_parameters(1);
        PredictedResponse.kappa = astrocyte_parameters(2);
        PredictedResponse.tau_d = astrocyte_parameters(3);     
        
        norm_fac = 1/PredictedResponse.tau_d;      % Keep area under the curve = 1
        H_t = (t>0).*(t/PredictedResponse.tau_d).*exp(-t/PredictedResponse.tau_d)*norm_fac;       % low-pass filter 1
        H_t(t<0) = 0;        % in case there are NaNs before t==0
        H_t(isnan(H_t)) = 0; 
    case 'low-pass_filter2'  % low-pass filter te^-t^2/2dt^2      
        PredictedResponse.w_f = astrocyte_parameters(1);
        PredictedResponse.kappa = astrocyte_parameters(2);
        PredictedResponse.tau_d = astrocyte_parameters(3);     
        
        norm_fac = 1/PredictedResponse.tau_d;      % Keep area under the curve = 1
        H_t = (t>0).*(t/PredictedResponse.tau_d).*exp(-t.^2/(2*PredictedResponse.tau_d^2))*norm_fac;       % low-pass filter 2
        H_t(t<0) = 0;        % in case there are NaNs before t==0
        H_t(isnan(H_t)) = 0;
                
    case 'low-pass_filter3'  % low-pass filter t^2e^-t^2/2dt^2      
        PredictedResponse.w_f = astrocyte_parameters(1);
        PredictedResponse.kappa = astrocyte_parameters(2);
        PredictedResponse.tau_d = astrocyte_parameters(3);     
        
        norm_fac = 4/(sqrt(pi)*PredictedResponse.tau_d);      % Keep area under the curve = 1
        H_t = (t>0).*(t.^2/PredictedResponse.tau_d^2).*exp(-t.^2/(PredictedResponse.tau_d^2))*norm_fac;       % low-pass filter 3
        H_t(t<0) = 0;        % in case there are NaNs before t==0
        H_t(isnan(H_t)) = 0;
end


%% Defining neural stimulus zeta(x,t) according to Aquino et al 2012 experiment
%  Gaussian in space and Block function in time

% FWHMx = 1;           % full width at half maximum of zeta(x)
% sigma_x = 1e-3*FWHMx/(2*sqrt(log(2)));       % standard deviation

sigma_x = 1e-3;             % standard deviation of Gaussian spread

zeta_x = exp(-xx.^2/(2*sigma_x^2));

stim_on_time = 8;           % how long the stimulus is on
stim_off_time = 12;      % how long the stimulus is off

zeta_t_vector = (t > 0).*(t <= stim_on_time) + ...
         (t > stim_on_time+stim_off_time).*(t <= 2*stim_on_time+stim_off_time) + ...
         (t > 2*stim_on_time+2*stim_off_time).*(t <= 3*stim_on_time+2*stim_off_time);
     
% Convolve zeta(t) with astrocytic function H_t to get final z(t)
% and make sure it has unit area
zeta_t_conv = conv(zeta_t_vector, H_t, 'same');
zeta_t_conv_norm = zeta_t_conv/trapz(t, zeta_t_conv);

% Transform zeta_t_conv into matrix 
zeta_t = repmat(zeta_t_conv_norm, length(x), 1);

% Get final z(x,t) = zeta_x * zeta_t_conv_norm
zeta_xt = zeta_x.*zeta_t;


%% Taking the Fourier transform of zeta(x,t)

zeta_kw = 1/p.Nk * ifftshift(fftshift(ifft(fft(zeta_xt, [], 1), [], 2), 1), 2);


%% Calculating the BOLD response 
% Y_xt is the total BOLD response
% t in x-axis and distance in y-axis

D  = p.rho_f*(2*Gamma - p.beta*p.Cz/p.tau);    
kz = sqrt((p.k_0)^2 + 1/(v_b)^2*p.Cz*(p.beta/p.tau)*(D/p.rho_f));

T_XiF =  p.Cz*p.rho_f*(D/p.rho_f - 1i*wmat)./(kmat.^2*v_b^2 + kz^2*v_b^2 - wmat.^2 - 2*1i*Gamma*wmat + eps);
T_FZeta = 1./(-(wmat + 0.5*1i*PredictedResponse.kappa).^2 + PredictedResponse.w_f^2);
T_QXi = (p.Q_0/p.Eps_0)*(1./(-1i*wmat + p.eta + p.tau^-1)).* ...
        (-p.V_0*1i*wmat + p.Cz*(p.eta - (p.tau^-1)*(p.beta - 2)));
    
T_XiZeta = T_XiF.*T_FZeta;

T_YZeta = ((p.k2 - p.k3)/p.rho_f)* ...
          (1 - (p.Eps_0/p.Q_0)*((p.k1 + p.k2)/(p.k2 - p.k3))*T_QXi).* ...
          T_XiZeta;

wvals = (-1).^(1:length(w));
kvals = (-1).^(1:length(k));

wvals_mat = repmat(wvals, length(k), 1);        %% prepare matrix for Fourier transform
kvals_mat = repmat(kvals.', 1, length(w));

Y_kw = T_YZeta.*zeta_kw;

% Note that x axis of matrix is time and y axis of matrix is distance
Y_xt = fftshift(kvals_mat.*ifft(kvals_mat.*wvals_mat.*fft(wvals_mat.*Y_kw,[],2),[],1));


%% If we want to response of the decomposition components
% Call T.Tx instead of T.Ttotal

% [~, ~, T, kz] = PoleDecomposition_num(p, v_b, Gamma, ...
%                 PredictedResponse.wf, PredictedResponse.kappa, kmat, wmat);
% 
% wvals = (-1).^(1:length(w));
% kvals = (-1).^(1:length(k));
% 
% wvals_mat = repmat(wvals, length(k), 1);        %% prepare matrix for Fourier transform
% kvals_mat = repmat(kvals.', 1, length(w));
% 
% Y_kw = T.Ttotal.*zeta_kw;
% Y_xt = fftshift(kvals_mat.*ifft(kvals_mat.*wvals_mat.*fft(wvals_mat.*Y_kw,[],2),[],1));

%% Just passing all resulting variables to PredictedSignal

% .' is used to invert matrix such that x-axis is distance and y-axis is time
PredictedResponse.Original_zeta_xt = (zeta_x.*repmat(zeta_t_conv, length(x), 1));
PredictedResponse.Norm_zeta_xt = zeta_xt;
PredictedResponse.Y_xt = Y_xt;

PredictedResponse.space = x;
PredictedResponse.time = t;



% 
% nt1 = find(t >= -5, 1 );
% %     nt2 = find(time >= 15, 1 );
% nt2 = find(t >= 25, 1 );
% 
% nx1 = find(x >= -0.01, 1 );
% nx2 = find(x >= 0.01, 1 );
% 
% Y_xt = Y_xt/(max(max(real(Y_xt(nx1:nx2,nt1:nt2)))));
% 
% figure;
% %         subplot(3,2,[1 3 5])
%         
% axes('fontsize',24','FontWeight','bold');
% 
% contourf(x(nx1:nx2)*1e3,t(nt1:nt2),real(Y_xt(nx1:nx2,nt1:nt2)).');grid on;%caxis([-0.4 0.32]);%colorbar;
% xlabel('x (mm)');ylabel('t (s)');
% colorbar;set(gca,'fontSize',24);
% xlim([-5,5]);
% ylim([-5,20]); 
% 
% figure;
% %         subplot(3,2,[1 3 5])
%         
% axes('fontsize',24','FontWeight','bold');
% contourf(x(nx1:nx2)*1e3,t(nt1:nt2),zeta_xt(nx1:nx2,nt1:nt2).');grid on;  % not normalized
% xlabel('x (mm)');ylabel('t (s)');
% colorbar;set(gca,'fontSize',24);
% xlim([-5,5]);
% ylim([0,15]); 
% ylim([0,20]);



    



