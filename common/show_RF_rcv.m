function show_RF_rcv(rf,label,fs,t0)

switch nargin
    case 1
        t = (0:size(rf,1)-1);
        label = 'RF plot';
    case 2
        t = (0:size(rf,1)-1);
    case 3
        t = (0:size(rf,1)-1)/fs;
    case 4
        t = (0:size(rf,1)-1)/fs+t0;
end
rf = rf./max(max(rf));
N_el = size(rf,2);

figure; hold on
for j = 1:N_el
    plot(t,rf(:,j)+j);
end
hold off
set(gca,'YDir','reverse')
title(label)
grid on
ylim([1 N_el])
ylabel('Rcv Channel #')
xlabel('Time')