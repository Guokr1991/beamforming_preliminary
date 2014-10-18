clear all; clc
addpath('/home/wjl11/Documents/MATLAB/field_ii');
addpath('/home/wjl11/Documents/MATLAB/field_ii_supplementary');
field_init(-1);

% define transducer
f0 = 5e6;
BW = 0.7;
N_el = 128;
elv_focus = 2;
elv_fnum = 6.5*elv_focus/0.04;
kerf_fraction = 0.05;
focus = [0 0 elv_focus];

% general parameters
c = 1540; fs = 100e6;

% derived parameters
lambda = c/f0;
el_height = elv_focus/elv_fnum;
el_width = (1-kerf_fraction)*lambda;
el_kerf = kerf_fraction*lambda;
n_sub_x = ceil(el_width/(lambda/2));
n_sub_y = ceil(el_height/(lambda/2));
el_pitch=el_width+el_kerf;
rx_pos = (-((N_el/2-1)*el_pitch+el_pitch/2):el_pitch:((N_el/2-1)*el_pitch+el_pitch/2));

% create transmit and receive arrays
Th=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);

% set impulse responses
tc=gauspuls('cutoff',f0,BW,-6,-40);
imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
xdc_impulse(Th,imp_resp);

position_pts = [0 0 0.06];
amplitude_pts = 100;

Ncyc = 1;
t_ex = 0:1/fs:Ncyc/f0;
pulse_ex = sin(2*pi*f0*t_ex);
xdc_excitation(Th, pulse_ex);
[tmp_rf, tmp_st] = calc_scat_multi(Th, Th, position_pts, amplitude_pts);
show_RF_rcv(tmp_rf);

% hold on
% plot(fliplr(tmp_rf(:,1))')
% plot(fliplr(tmp_rf(:,50)'),'r')
% hold off
% 
% f0=3e6; % Transducer center frequency [Hz]
% fs=100e6; % Sampling frequency [Hz]
% c=1540; % Speed of sound [m/s]
% lambda=c/f0; % Wavelength [m]
% height=5/1000; % Height of element [m]
% width=1/1000; % Width of element [m]
% kerf=width/4; % Distance between transducer elements [m]
% N_elements=32; % Number of elements
% focus=[0 0 60]/1000; % Initial electronic focus
% % Define the transducer
% Th = xdc_linear_array (N_elements, width, height, kerf, 2, 3, focus);
% % Set the impulse response and excitation of the emit aperture
% 
% impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
% impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
% xdc_impulse (Th, impulse_response);
% excitation=sin(2*pi*f0*(0:1/fs:2/f0));
% xdc_excitation (Th, excitation);
% % Do the calculation
% [v,t]=calc_scat_multi (Th, Th, [0 0 40]/1000, 1);
% show_RF_rcv(v);
% % Plot the individual responses
% subplot(211)
% [N,M]=size(v);
% v=v/max(max(v));
% for i=1:N_elements
% plot((0:N-1)/fs+t,v(:,i)+i), hold on
% end
% hold off
% title('Individual traces')
% xlabel('Time [s]')
% ylabel('Normalized response')
% subplot(212)
% plot((0:N-1)/fs+t,sum(v'))
% title('Summed response')
% xlabel('Time [s]')
% ylabel('Normalized response')