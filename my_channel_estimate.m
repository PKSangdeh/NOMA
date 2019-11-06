%% My channel Estimation for DL of M single-antenna users
% The single-antenna AP sends a signal, 
% M users receives the frame and estimate the channel coefficient
%%
clc
close all
loc=1;
show_detail=0;
%% Parameters
num_users = 3;
num_valid_sc = 52;
num_subcarrier = 64;
num_symbol_per_frame = 20;
num_frame = 20;
num_sample_process = 50000;
num_sample_shift = 2e6;
num_move_sample = 0;
M=2;
%% Preamble info
% STF
STF_POS = [3  7  11  15  19  23  30  34  38  42  46  50];
STF_DAT_GRP = zeros(1, 52);    
STF_DAT_GRP(1,STF_POS) = sqrt(2)*[-1-1i  -1+1i  -1-1i  -1-1i  -1-1i  1-1i  -1-1i  -1-1i  1+1i  1-1i  1+1i  -1+1i];
% LTF
LTF_DAT_GRP = zeros(2, 52);
LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(2,:) = [1 1 -1 -1 1 1 -1 -1 1 1 -1 -1 -1 1 -1 1 1 1 -1 -1 -1 -1 1 1 1 -1 1 1 1 -1 -1 -1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 -1 -1 1 -1 -1 1 1 -1];
% Pilot
POS_VALID_SC = [39:64 2:27];
POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
POS_PILOT_SC = [6 20 33 47];
PILOT_DAT_GRP = zeros(4,num_symbol_per_frame);
PILOT_DAT_GRP(1,:) = [-1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1];
PILOT_DAT_GRP(2,:) = [1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 1 1 -1 -1 1 1 1];
PILOT_DAT_GRP(3,:) = [-1 -1 1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 1 1 1];
PILOT_DAT_GRP(4,:) = [1 1 1 -1 1 1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 1 1 -1];
% Preamble
PREAMABLE = [1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,1,1,-1,1,-1,-1,1,-1,1];
PREAMABLE = [PREAMABLE PREAMABLE];
%% Frame detection, synchornization, and channel estimation
%% read received signals from files
signal_stream_original = zeros(num_users, num_sample_process);
for ii=1:num_users
fileID = fopen(['C:\Users\pedra\Desktop\NOMA_2X3\signals\signal_rx_u' num2str(ii) '.dat'], 'rb');
    if (fileID < 0)
        disp('Error: fail to open files!');
        pause;
    end
    frewind(fileID);
    fseek(fileID, num_sample_shift, 'bof');
    data_ant_f = fread(fileID, 2*num_sample_process, 'float');
    data_ant_c = transpose(data_ant_f(1:2:end) + 1i*data_ant_f(2:2:end));
    signal_stream_original(ii,:) = data_ant_c(1:num_sample_process);
    fclose(fileID);
end
signal_stream_combined = signal_stream_original;
%% time synchornization
cross_correlation_arr = zeros(num_users,num_sample_process);
data_freq = zeros(1, num_subcarrier);
data_freq(POS_VALID_SC) = LTF_DAT_GRP(1,:);
data_time = ifft(data_freq)*sqrt(num_subcarrier);
ltf_signal = [data_time(end-15:end) data_time];
local_signal = [ltf_signal  ltf_signal];
for ii = 1:num_sample_process-length(local_signal)+1
    for jj=1:num_users
    signal_segment=signal_stream_combined(jj,ii:ii+length(local_signal)-1);
    cross_correlation_arr(jj,ii) =  local_signal*signal_segment'/(norm(local_signal)*norm(signal_segment));
    end
end
if show_detail~=0
figure; hold on;
for i=1:num_users
    subplot(num_users,1,i)
    plot(abs(cross_correlation_arr(i,:)),'b-');
    xlabel('sample index');
    ylabel('correlation');
    title('Cross correlation (DL)');
    box on;grid on; 
    axis([0 length(signal_stream_combined) 0 1]);
end
end
n_symb=80*num_symbol_per_frame;
chan_coeff_matrix = zeros(num_users,M,num_valid_sc); 
for u_idx=1:num_users
    %% frame capture
    [~, begofframe] = max(abs(cross_correlation_arr(u_idx,1:1600)));
    begofframe = begofframe + n_symb - 160 - num_move_sample;
    captured_sig=signal_stream_combined(u_idx,begofframe:end);
    %% freq synch
    freq_offset1 = 0;
    freq_offset2 = 0;
    for ii = 0:1600:length(captured_sig)-1600
        % CP_based
        for jj = 1:20
            freq_offset1 = freq_offset1 + captured_sig(ii+80*(jj-1)+1:ii+80*(jj-1)+16)*captured_sig(ii+80*(jj-1)+65:ii+80*(jj-1)+80)';
        end
        % LTF_based
        freq_offset2 = freq_offset2 + captured_sig(ii+177:ii+240)*captured_sig(ii+257:ii+320)';
    end
    angle_os1 = angle(freq_offset1)/64;
    angle_os2 = angle(freq_offset2)/80;
    freq_comp = exp(1i*angle_os1*(1:length(captured_sig)));
    freq_comp_fram = repmat(freq_comp, size(captured_sig,1), 1);
    signal_stream_buf = captured_sig.*freq_comp_fram;
    %disp(['freq offset User ' num2str(u_idx) '-  CP_based:' num2str(angle_os1)  '  LTF_based:' num2str(angle_os2)]);
    %% Channel estimation
        % capture one frame
        frm_idx=1;
        signal_frame=signal_stream_buf(1600*(frm_idx-1)+1:1600*frm_idx);
        data_t1 = signal_frame;
        data_t2 = reshape(data_t1, num_subcarrier+16, num_symbol_per_frame);
        data_t3 = data_t2(17:end,:);
        data_f1 = fft(data_t3)/sqrt(num_subcarrier);
        signal_matrix = data_f1(POS_VALID_SC,:);
        % estimate channel
        signal_ltf_avg = (signal_matrix(:, 3) + signal_matrix(:, 4))/2;
        chan_coeff_matrix(u_idx,1,:) = signal_ltf_avg./LTF_DAT_GRP(1, :)';
        signal_ltf_avg2 = (signal_matrix(:, 5) + signal_matrix(:, 6))/2;
        chan_coeff_matrix(u_idx,2,:) = signal_ltf_avg2./LTF_DAT_GRP(2, :)';
end
%% Plot Channel Responses
for u_idx = 1:num_users
    figure; 
    subplot(2, 2, 1);
    title(['Magnitude of h_1_' num2str(u_idx)]);    
    hold on;
    grid on;    
    col=abs(rand(1,3));
    plot(abs(squeeze(chan_coeff_matrix(u_idx,1,:))),'*-','Color',col);
    axis([1 num_valid_sc 0 2*max(abs(chan_coeff_matrix(u_idx,1,:)))])
    
    subplot(2, 2, 2);
    title(['Phase of h_1_' num2str(u_idx)]);    
    hold on;
    grid on;    
    plot(unwrap(angle((squeeze(chan_coeff_matrix(u_idx,1,:))))),'*-','Color',col);
    
    subplot(2, 2, 3);
    title(['Magnitude of h_2_' num2str(u_idx)]);    
    hold on;
    grid on;    
    plot(abs(squeeze(chan_coeff_matrix(u_idx,2,:))),'*-','Color',col);
    axis([1 num_valid_sc 0 2*max(abs(chan_coeff_matrix(u_idx,2,:)))])
    
    subplot(2, 2, 4);
    title(['Phase of h_2_' num2str(u_idx)]);    
    hold on;
    grid on;    
    plot(unwrap(angle(squeeze((chan_coeff_matrix(u_idx,2,:))))),'*-','Color',col)
end
ch_u1=squeeze(chan_coeff_matrix(1,:,:));
ch_u2=squeeze(chan_coeff_matrix(2,:,:));
ch_u3=squeeze(chan_coeff_matrix(3,:,:));
save('C:\Users\pedra\Desktop\NOMA_2X3\channel_u1.mat','ch_u1'); 
save('C:\Users\pedra\Desktop\NOMA_2X3\channel_u2.mat','ch_u2');
save('C:\Users\pedra\Desktop\NOMA_2X3\channel_u3.mat','ch_u3');
