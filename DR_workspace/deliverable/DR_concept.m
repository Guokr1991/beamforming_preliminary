function [focus_rf, x, z] = DR_concept(rf,acq_params,bf_params)
% Dynamic receive code - Will Long. Latest revision: 10/1/14
% Inputs: 
% rf - raw rf data organized [rf_line,rx_chan, tx_event]
% acq_params - parameters include rx_pos, c, t0
% bf_params - parameters include x (tx_pos or A line lateral location)
% 
% NOTE: for the case when all rx positions are relative and identical for
% each tx location (i.e. image using 128 subaperture of 256 array with tx
% focus at center of subaperture. Walk subaperture to generate multiple
% scan lines. In this case, lateral distance from focus is only the rx pos.


%all focusing performed on axis

% acq_params.rx_pos; %lateral positions of each rx element relative to tx
% bf_params.x; %lateral positions of each tx focus (per tx event)
% acq_params.c;
% acq_params.t0; start time of rf data reference in terms of sample number
tic

x = bf_params.x;
% z_ref_test = single(acq_params.t0+1:size(rf,1)/acq_params.fs)*acq_params.c;
z_ref = ((acq_params.t0+1:acq_params.t0+size(rf,1))/acq_params.fs)*acq_params.c;
z = z_ref/2;

[b,a]=butter(2,[.05 .95]);

% intialize matrices and arrays for speed
n_tx = length(bf_params.x);
n_rcv_chn = length(acq_params.rx_pos);
n_depth = length(z);
focus_rf = zeros(length(z),n_tx);

if n_rcv_chn ~= size(rf,2) || n_tx ~= size(rf,3)
    disp('Mismatch in RF data.')
end
    
dz = repmat(z',1,n_rcv_chn);
dx = repmat(acq_params.rx_pos,n_depth,1);
dr = sqrt(dz.^2+dx.^2);
t_samp = (dr+dz)./acq_params.c;

for j = 1:n_tx % iterate for every tx event
    focus_rf(:,j) = sum(lin_interp(z_ref'/acq_params.c,squeeze(rf(:,:,j)),t_samp),2);
end

focus_rf = filter(b,a,focus_rf);
focus_rf(find(isnan(focus_rf))) = 0; 
toc

