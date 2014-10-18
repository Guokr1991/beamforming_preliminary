% rf - raw rf data organized [rf_line,rx_chan, tx_event]
function [focus_rf, x, z] = dynamic_receive_concept(rf,acq_params,bf_params)

%all focusing performed on axis

% acq_params.rx_pos; %lateral positions of each rx element relative to tx
% bf_params.x; %lateral positions of each tx focus (per tx event)
% acq_params.c;
% acq_params.t0;

x = bf_params.x;
z_ref = single((acq_params.t0+1:acq_params.t0+size(rf,1))/acq_params.fs)...
    *acq_params.c;
z = z_ref/2;

% intialize matrices and arrays for speed
n_tx = length(bf_params.x);
n_rcv_chn = length(acq_params.rx_pos);
focus_rf = zeros(length(z),n_tx);
t_samp = zeros(1,n_rcv_chn);

if n_rcv_chn ~= size(rf,2) || n_tx ~= size(rf,3)
    disp('Mismatch in RF data.')
end

for j = 1:n_tx % iterate for every tx event
    for i = 1:length(z_ref) % iterate for each depth
        d_tx = z_ref(i);

        for k = 1:n_rcv_chn % iterate for each rx element
            d_rx = sqrt((acq_params.rx_pos(k)-bf_params.x(j))^2+d_tx^2);
            t_samp(k) = (d_tx+d_rx)/acq_params.c;
        end

        focus_rf(i,j) = sum(sum(eye(n_rcv_chn).*...
            interp1(z_ref'/acq_params.c,squeeze(rf(:,:,j)),t_samp,...
            'linear'))); 
    end
    disp(['tx line ' num2str(j) ' processed.'])
end

