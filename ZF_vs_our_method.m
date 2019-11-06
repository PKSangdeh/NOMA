%% NOMA Stations
clear;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load alph
M=4;
num_ant_ap = 1;
num_stm_ap = 1;
num_sample_processing = 16000;
num_sc = 64;
num_symbol_per_frame = 20;
num_sample_shift = 1e6;
num_averaging_sc = 4;
num_valid_sc = 52;
const_num_sc = 64;
const_correl_threhold = 0.5;
subc_smooth_enabled=0;
minimal=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edg_ign=0;
show_detail=1; % shows every details and steps
re_synch_enable=1; % freq. re-synchornization after subtraction
perv_noma=1;       % Our previous version of NOMA
bypass_freq_synch=0;
bypass_phase_corr=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STF 
STF_POS = [3  7  11  15  19  23  30  34  38  42  46  50];
STF_DAT_GRP = zeros(2, 52);    
STF_DAT_GRP(1,STF_POS) = sqrt(2)*[-1-1i  -1+1i  -1-1i  -1-1i  -1-1i  1-1i  -1-1i  -1-1i  1+1i  1-1i  1+1i  -1+1i];
STF_DAT_GRP(2,STF_POS) = sqrt(2)*[-1+1i  -1+1i  1-1i  1-1i  -1-1i  1-1i  -1+1i  1+1i  1-1i  -1+1i  1-1i  -1+1i]; 
%%% LTF
LTF_DAT_GRP = zeros(2, 52);
LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(2,:) = [1 1 -1 -1 1 1 -1 -1 1 1 -1 -1 -1 1 -1 1 1 1 -1 -1 -1 -1 1 1 1 -1 1 1 1 -1 -1 -1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 -1 -1 1 -1 -1 1 1 -1];
%%%% pilot
POS_VALID_SC = [39:64 2:27];
POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
POS_PILOT_SC = [6 20 33 47];
PILOT_DAT_GRP = zeros(4,num_symbol_per_frame);
PILOT_DAT_GRP(1,:) = [-1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1];
PILOT_DAT_GRP(2,:) = [1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 1 1 -1 -1 1 1 1];
PILOT_DAT_GRP(3,:) = [-1 -1 1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 1 1 1];
PILOT_DAT_GRP(4,:) = [1 1 1 -1 1 1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 1 1 -1];
% PILOT_DAT_GRP = 1/sqrt(2)*ones(4,num_symbol_per_frame)+1/sqrt(2)*ones(4,num_symbol_per_frame)*1j;

%% weak user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_id = fopen('signals/noma_rx_u2.dat', 'rb');%
if (file_id < 0)
    error('Error: fail to open files!');
end
frewind(file_id);
fseek(file_id, num_sample_shift, 'bof');
rx_signal_f = fread(file_id, 2*num_sample_processing, 'float');
rx_signal_buf = transpose(rx_signal_f(1:2:end) + 1i*rx_signal_f(2:2:end));
fclose(file_id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Cross correlation
ltf_freq_signal = zeros(1, num_sc);
ltf_freq_signal(POS_VALID_SC) = (1-alph)*LTF_DAT_GRP(1,:)+alph*LTF_DAT_GRP(2,:);
lte_time_signal = ifft(ltf_freq_signal)*sqrt(num_sc);
ltf_waveform = [lte_time_signal(end-15:end) lte_time_signal lte_time_signal(end-15:end) lte_time_signal];
cross_correl_value = zeros(1,num_sample_processing);
for ii = 161:num_sample_processing-160
    signal_segment = rx_signal_buf(ii:ii+160-1);
    cross_correl_value(ii) =  ltf_waveform*signal_segment'/(norm(ltf_waveform)*norm(signal_segment));
end
[max_correl_value, max_correl_pos] = max(abs(cross_correl_value(1:end-1600)));
begofframe = max_correl_pos - 160;
if show_detail==1 && minimal==0
   figure(); hold on; grid on; box on;
   plot(abs(cross_correl_value),'k-');
   xlabel('sample index');
   ylabel('correlation');
   title('Cross-correlation on weak user')
   axis([0 length(rx_signal_buf) 0 1]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the begining of a frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_signal_buf = rx_signal_buf(begofframe:end);
rx_signal_frame = rx_signal_buf(1:1600);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_offset_cp_based = 0;
freq_offset_ltf_based = 0;
% CP
for jj = 1:20
    freq_offset_cp_based = freq_offset_cp_based + rx_signal_frame(80*(jj-1)+1:80*(jj-1)+16)*rx_signal_frame(80*(jj-1)+65:80*(jj-1)+80)';
end
% LTF
freq_offset_ltf_based = freq_offset_ltf_based + rx_signal_frame(177:240)*rx_signal_frame(257:320)';
delta_phase_per_sample_cp_based = angle(freq_offset_cp_based)/64;
delta_phase_per_sample_ltf_based = angle(freq_offset_ltf_based)/80;
delta_phase_per_sample = delta_phase_per_sample_ltf_based;
rx_dds = exp(1i*delta_phase_per_sample*(1:size(rx_signal_frame,2)));
rx_dds_frame = repmat(rx_dds, size(rx_signal_frame,1), 1);
rx_signal_frame = rx_signal_frame.*rx_dds_frame;
if show_detail==1 
   disp('Weak User')
   disp(['cp-based freq offset estimate: ' num2str(delta_phase_per_sample_cp_based)]);
   disp(['ltf-based freq offset estimate: ' num2str(delta_phase_per_sample_ltf_based)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% received frequency frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_signal_one_frame = rx_signal_frame;
freq_signal_one_frame = zeros(num_ant_ap, num_valid_sc, num_symbol_per_frame);
for rx_idx = 1:num_ant_ap
    time_signal_tmp1 = time_signal_one_frame(rx_idx,:);
    time_signal_tmp2 = reshape(time_signal_tmp1, num_sc+16, num_symbol_per_frame);
    time_signal_tmp3 = time_signal_tmp2(17:end,:);
    freq_signal_tmp1 = fft(time_signal_tmp3)/sqrt(num_sc);
    freq_signal_tmp2 = freq_signal_tmp1(POS_VALID_SC,:);
    freq_signal_one_frame(rx_idx,:,:) = freq_signal_tmp2;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Intrf resilient Receiver
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute G
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute EYY and EYX
    EYY_SUM = zeros(num_ant_ap, num_ant_ap, num_valid_sc);
    EYX_SUM = zeros(num_ant_ap, num_valid_sc);
    ltf_dat=(1-alph)*LTF_DAT_GRP(1,:)+alph*LTF_DAT_GRP(2,:);
    stf_dat=(1-alph)*STF_DAT_GRP(1,:)+alph*STF_DAT_GRP(2,:);
    freq_frame=squeeze(freq_signal_one_frame);
    for sc_idx = 1:num_valid_sc
        % TSF (OFDM 1)
        Y = freq_signal_one_frame(:, sc_idx, 1);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = stf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % TSF (OFDM 2)
        Y = freq_signal_one_frame(:, sc_idx, 2);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = stf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % LTF (OFDM 3)
        Y = freq_signal_one_frame(:, sc_idx, 3);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = ltf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % LTF (OFDM 4)
        Y = freq_signal_one_frame(:, sc_idx, 4);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = ltf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
    end
    % compute G
    G_ARR_PER_FRAME = zeros(num_ant_ap, num_valid_sc);
    for sc_idx = 1:num_valid_sc
        EYY = zeros(num_ant_ap, num_ant_ap);
        EYX = zeros(num_ant_ap, 1);
        sc_lwbd = max(1, sc_idx-num_averaging_sc);
        sc_upbd = min(sc_idx+num_averaging_sc, num_valid_sc);
        for kk = sc_lwbd:sc_upbd
            EYY = EYY + squeeze(EYY_SUM(:,:,kk));
            EYX = EYX + squeeze(EYX_SUM(:,kk));
        end
        G = pinv(EYY)*(EYX);
        G_ARR_PER_FRAME(:,sc_idx) = G;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decoded data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_data_frame = zeros(num_valid_sc,num_symbol_per_frame);
    for sc_idx = 1:num_valid_sc
        G = G_ARR_PER_FRAME(:, sc_idx);
        Y = squeeze(freq_signal_one_frame(:, sc_idx, :));
        decoded_data_frame(sc_idx,:) = G'*Y;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % phase calibration for decoded data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phase_offset1 = decoded_data_frame(POS_PILOT_SC(1),:)./PILOT_DAT_GRP(1,:);
    phase_offset2 = decoded_data_frame(POS_PILOT_SC(2),:)./PILOT_DAT_GRP(2,:);
    phase_offset3 = decoded_data_frame(POS_PILOT_SC(3),:)./PILOT_DAT_GRP(3,:);
    phase_offset4 = decoded_data_frame(POS_PILOT_SC(4),:)./PILOT_DAT_GRP(4,:);
    decoded_data_frame(1:13,:) = decoded_data_frame(1:13,:)./repmat(phase_offset1,13,1);
    decoded_data_frame(14:27,:) = decoded_data_frame(14:27,:)./repmat(phase_offset2,14,1);
    decoded_data_frame(28:40,:) = decoded_data_frame(28:40,:)./repmat(phase_offset3,13,1);
    decoded_data_frame(41:52,:) = decoded_data_frame(41:52,:)./repmat(phase_offset4,12,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decoded data with phase compension
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_payload_w_phase_comp = decoded_data_frame(POS_PAYLOAD_SC,5:num_symbol_per_frame);  
decoded_payload_w_phase_comp_arr = decoded_payload_w_phase_comp(:);       
rx_data=decoded_payload_w_phase_comp_arr;
demoded_data=qamdemod(rx_data,M,'UnitAveragePower',true); 
ref_points=qammod(demoded_data,M,'UnitAveragePower',true); 
err=abs(rx_data-ref_points);
evm=10*log10(mean(err.^2));
h = figure(); 
hold on; grid on;box on;
title(['constellation at weak user, EVM=',num2str(evm)]);            
        scatter(real(decoded_payload_w_phase_comp_arr),imag(decoded_payload_w_phase_comp_arr),'ko','linewidth',0.5,'MarkerFaceColor','b');
        axis([-2, 2, -2, 2]);
          
        
 %%       
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % %%%%%%%%%%%%%%%% STRONG USER %%%%%%%%%%%%%%% STRONG USER %%%%%%%%%%% STRONG USER %%%%%%%%%%%
 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_id = fopen('signals/noma_rx_u1.dat', 'rb');%
if (file_id < 0)
   error('Error: fail to open files!');
end
frewind(file_id);
fseek(file_id, num_sample_shift, 'bof');
rx_signal_f = fread(file_id, 2*num_sample_processing, 'float');
rx_signal_buf = transpose(rx_signal_f(1:2:end) + 1i*rx_signal_f(2:2:end));
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Cross correlation
ltf_freq_signal = zeros(1, num_sc);
ltf_freq_signal(POS_VALID_SC) =(1-alph)*LTF_DAT_GRP(1,:)+alph*LTF_DAT_GRP(2,:);
lte_time_signal = ifft(ltf_freq_signal)*sqrt(num_sc);
ltf_waveform = [lte_time_signal(end-15:end) lte_time_signal lte_time_signal(end-15:end) lte_time_signal];
cross_correl_value = zeros(1,num_sample_processing);
for ii = 161:num_sample_processing-160
    signal_segment = rx_signal_buf(ii:ii+160-1);
    cross_correl_value(ii) =  ltf_waveform*signal_segment'/(norm(ltf_waveform)*norm(signal_segment));
end
[max_correl_value, max_correl_pos] = max(abs(cross_correl_value(1:end-1600)));
begofframe = max_correl_pos - 160;
if show_detail==1 && minimal==0
   figure(); hold on; grid on; box on;
   plot(abs(cross_correl_value),'k-');
   xlabel('sample index');
   ylabel('correlation');
   axis([0 length(rx_signal_buf) 0 1]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the begining of a frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_signal_buf = rx_signal_buf(begofframe:end);
rx_signal_frame = rx_signal_buf(1:1600);
rx_original=reshape(rx_signal_frame,80,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_offset_cp_based = 0;
freq_offset_ltf_based = 0;
% CP
for jj = 1:20
    freq_offset_cp_based = freq_offset_cp_based + rx_signal_frame(80*(jj-1)+1:80*(jj-1)+16)*rx_signal_frame(80*(jj-1)+65:80*(jj-1)+80)';
end
% LTF
freq_offset_ltf_based = freq_offset_ltf_based + rx_signal_frame(177:240)*rx_signal_frame(257:320)';
delta_phase_per_sample_cp_based = angle(freq_offset_cp_based)/64;
delta_phase_per_sample_ltf_based = angle(freq_offset_ltf_based)/80;
delta_phase_per_sample = delta_phase_per_sample_ltf_based;
rx_dds = exp(1i*delta_phase_per_sample*(1:size(rx_signal_frame,2)));
if bypass_freq_synch==0
    rx_dds_frame = repmat(rx_dds, size(rx_signal_frame,1), 1);
else
    rx_dds_frame = ones(1,1600);
    disp('1st freq. synchornization is bypassed')
end
rx_signal_frame = rx_signal_frame.*rx_dds_frame;
if show_detail==1 && bypass_freq_synch==0
   fprintf('\n')
   disp('Strong User')
   disp(['cp-based freq offset estimate: ' num2str(delta_phase_per_sample_cp_based)]);
   disp(['ltf-based freq offset estimate: ' num2str(delta_phase_per_sample_ltf_based)]);
end


%>>>>>>>>>>>>> ref signal
% ref signal
ref_signal=reshape(rx_signal_frame,80,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% received frequency frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_signal_one_frame = rx_signal_frame;
freq_signal_one_frame = zeros(num_ant_ap, num_valid_sc, num_symbol_per_frame);
for rx_idx = 1:num_ant_ap
    time_signal_tmp1 = time_signal_one_frame(rx_idx,:);
    time_signal_tmp2 = reshape(time_signal_tmp1, num_sc+16, num_symbol_per_frame);
    time_signal_tmp3 = time_signal_tmp2(17:end,:);
    freq_signal_tmp1 = fft(time_signal_tmp3)/sqrt(num_sc);
    freq_signal_tmp2 = freq_signal_tmp1(POS_VALID_SC,:);
    freq_signal_one_frame(rx_idx,:,:) = freq_signal_tmp2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_t1 = rx_signal_frame;
data_t2 = reshape(data_t1, num_sc+16, []);
data_t3 = data_t2(17:end,:);
data_f1 = fft(data_t3)/sqrt(num_sc);
signal_matrix = data_f1(POS_VALID_SC,:);
% estimate channel
signal_ltf_avg = (signal_matrix(:, 3) + signal_matrix(:, 4))/2;
ltf_dat=(1-alph)*LTF_DAT_GRP(1,:)+alph*LTF_DAT_GRP(2,:);
chan_coeff_matrix=signal_ltf_avg./ltf_dat';
if show_detail==1 && minimal==0
    figure
    plot(1:num_valid_sc,abs(chan_coeff_matrix),'-ok','linewidth',2,'MarkerFaceColor','r')
    xlabel('Subcarrier index')
    ylabel('absolute magnitude')
    title('Estimated channel at strong user')
    axis([1 52 0 1.5*max(abs(chan_coeff_matrix))])
    grid on
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>> ZF-NOMA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       
% ZF         
% Subcarrier smoothing
ch_smooth=zeros(52,1);
ch=[chan_coeff_matrix(1);chan_coeff_matrix;chan_coeff_matrix(52)];
for ii=2:53
ch_smooth(ii-1)=(ch(ii)+ch(ii-1)+ch(ii+1))/3;
end
if subc_smooth_enabled==1
chan_coeff_matrix=ch_smooth;
end
% Received Frame with no phase estimation
decoded_data_frame_zf = squeeze((freq_signal_one_frame))./repmat(chan_coeff_matrix,1,20);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % phase calibration for decoded data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_offset1 = decoded_data_frame_zf(POS_PILOT_SC(1),5:end)./PILOT_DAT_GRP(1,5:end);
phase_offset2 = decoded_data_frame_zf(POS_PILOT_SC(2),5:end)./PILOT_DAT_GRP(2,5:end);
phase_offset3 = decoded_data_frame_zf(POS_PILOT_SC(3),5:end)./PILOT_DAT_GRP(3,5:end);
phase_offset4 = decoded_data_frame_zf(POS_PILOT_SC(4),5:end)./PILOT_DAT_GRP(4,5:end);

phase_offset=ones(52,20);
phase_offset(1:13,5:end) =repmat(phase_offset1,13,1);
phase_offset(14:27,5:end)=repmat(phase_offset2,14,1);
phase_offset(28:40,5:end)=repmat(phase_offset3,13,1);
phase_offset(41:52,5:end)=repmat(phase_offset4,12,1);
phase_corr=1./phase_offset;
decoded_data_frame_zf=decoded_data_frame_zf.*phase_corr;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decoded data with phase compension
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_payload_w_phase_comp_zf = decoded_data_frame_zf(POS_PAYLOAD_SC,5:num_symbol_per_frame);  
    decoded_payload_w_phase_comp_arr_zf = decoded_payload_w_phase_comp_zf(:);
    rx_data_zf=decoded_payload_w_phase_comp_arr_zf;
    ref_points=qammod(qamdemod(rx_data_zf,M,'UnitAveragePower',true),M,'UnitAveragePower',true);
    err=abs(rx_data_zf-ref_points);
    evm=10*log10(mean(err.^2));
    figure();
    hold on; grid on;box on;
    title(['constellation at Strong user-ZF, EVM=',num2str(evm)]);
    hold on; grid on;box on;          
    scatter(real(decoded_payload_w_phase_comp_arr_zf),imag(decoded_payload_w_phase_comp_arr_zf),'ko','linewidth',0.5,'MarkerFaceColor',[0.1 0.3 0.1]);
    axis([-2, 2, -2, 2]);
    
    rx_data_p=decoded_payload_w_phase_comp_arr_zf;   
    demoded_interf=qamdemod(rx_data_p,M,'UnitAveragePower',true);   
    interf_p=qammod(demoded_interf,M,'UnitAveragePower',true);
    desired_sig=rx_data_p-interf_p*alph;
    desired_sig=desired_sig/(1-alph);
    demoded_data=qamdemod(desired_sig,M,'UnitAveragePower',true); 
    ref_points=qammod(demoded_data,M,'UnitAveragePower',true); 
    err=abs(ref_points-desired_sig);
    evm=10*log10(mean(err.^2));
    figure()
    scatter(real(desired_sig),imag(desired_sig),'ko','linewidth',0.5,'MarkerFaceColor',[0.1 0.2 0.1]);
    axis([ -2 2 -2 2])
    title(['Desired signal-ZF, EVM=',num2str(evm)])
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute EYY and EYX
EYY_SUM = zeros(num_ant_ap, num_ant_ap, num_valid_sc);
EYX_SUM = zeros(num_ant_ap, num_valid_sc);
ltf_dat=(1-alph)*LTF_DAT_GRP(1,:)+alph*LTF_DAT_GRP(2,:);
stf_dat=(1-alph)*STF_DAT_GRP(1,:)+alph*STF_DAT_GRP(2,:);
for sc_idx = 1:num_valid_sc
    % TSF (OFDM 1)
        Y = freq_signal_one_frame(:, sc_idx, 1);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = stf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % TSF (OFDM 2)
        Y = freq_signal_one_frame(:, sc_idx, 2);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = stf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % LTF (OFDM 3)
        Y = freq_signal_one_frame(:, sc_idx, 3);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = ltf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % LTF (OFDM 4)
        Y = freq_signal_one_frame(:, sc_idx, 4);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = ltf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
end
    % compute G
    G_ARR_PER_FRAME = zeros(num_ant_ap, num_valid_sc);
for sc_idx = 1:num_valid_sc
        EYY = zeros(num_ant_ap, num_ant_ap);
        EYX = zeros(num_ant_ap, 1);
        sc_lwbd = max(1, sc_idx-num_averaging_sc);
        sc_upbd = min(sc_idx+num_averaging_sc, num_valid_sc);
        for kk = sc_lwbd:sc_upbd
            EYY = EYY + squeeze(EYY_SUM(:,:,kk));
            EYX = EYX + squeeze(EYX_SUM(:,kk));
        end
        G = pinv(EYY)*(EYX);
        G_ARR_PER_FRAME(:,sc_idx) = G;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decoded data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
decoded_data_frame = zeros(num_valid_sc,num_symbol_per_frame);
for sc_idx = 1:num_valid_sc
        G = G_ARR_PER_FRAME(:, sc_idx);
        Y = squeeze(freq_signal_one_frame(:, sc_idx, :));
        decoded_data_frame(sc_idx,:) = G'*Y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase calibration for decoded data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bypass_freq_synch==0
phase_offset1 = decoded_data_frame(POS_PILOT_SC(1),5:end)./PILOT_DAT_GRP(1,5:end);
phase_offset2 = decoded_data_frame(POS_PILOT_SC(2),5:end)./PILOT_DAT_GRP(2,5:end);
phase_offset3 = decoded_data_frame(POS_PILOT_SC(3),5:end)./PILOT_DAT_GRP(3,5:end);
phase_offset4 = decoded_data_frame(POS_PILOT_SC(4),5:end)./PILOT_DAT_GRP(4,5:end);

phase_offset=ones(52,20);
phase_offset(1:13,5:end) =repmat(phase_offset1,13,1);
phase_offset(14:27,5:end)=repmat(phase_offset2,14,1);
phase_offset(28:40,5:end)=repmat(phase_offset3,13,1);
phase_offset(41:52,5:end)=repmat(phase_offset4,12,1);
phase_corr=1./phase_offset;
decoded_data_frame=decoded_data_frame.*phase_corr;
else
    phase_corr=ones(52,20);
    disp('1st Phase correction is bypassed')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decoded data with phase compension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
decoded_payload_w_phase_comp = decoded_data_frame(POS_PAYLOAD_SC,5:num_symbol_per_frame);
decoded_payload_w_phase_comp_arr = decoded_payload_w_phase_comp(:);
rx_data_1st=decoded_payload_w_phase_comp_arr;
ref_points=qammod(qamdemod(rx_data_1st,M,'UnitAveragePower',true),M,'UnitAveragePower',true);
err=abs(rx_data_1st-ref_points);
evm=10*log10(mean(err.^2));
figure();
hold on; grid on;box on;
title(['constellation at Strong user-G filter, EVM=',num2str(evm)]);
scatter(real(decoded_payload_w_phase_comp_arr),imag(decoded_payload_w_phase_comp_arr),'ko','linewidth',0.5,'MarkerFaceColor','r');
axis([-2, 2, -2, 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interference Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_data=decoded_payload_w_phase_comp_arr;
demoded_interf=qamdemod(rx_data,M,'UnitAveragePower',true);
interf=reshape(qammod(demoded_interf,M,'UnitAveragePower',true),48,[]);
stf_dat = STF_DAT_GRP(2,:);
ltf_dat = LTF_DAT_GRP(2,:);
data_freq = zeros(num_valid_sc,num_symbol_per_frame);
data_freq(:,1) = stf_dat.';
data_freq(:,2) = stf_dat.';
data_freq(:,3) = ltf_dat.';
data_freq(:,4) = ltf_dat.';
for jj = 5:num_symbol_per_frame
    data_freq(POS_PAYLOAD_SC,jj) = interf(:,jj-4);
    data_freq(POS_PILOT_SC,jj) = PILOT_DAT_GRP(:,jj);
end
data_f_frm2 = zeros(num_sc, num_symbol_per_frame);
data_f_frm2(POS_VALID_SC,:) = data_freq;
f_tx2=data_f_frm2;
intf=ifft(f_tx2)*sqrt(num_sc);
intf=[intf(end-15:end,:); intf];
load time_signal_tx2
% Verifying re-constructed interference
if sum(sum(intf-time_signal_tx2))==0
   fprintf('\n')
   disp('Reconstruction is verified')
   fprintf('\n')
else
   fprintf('\n')
   disp('There is at least one mismatch between real and reconstructed interference')
end
fprintf('Performance:\nZF approach : %f\n',evm)
chan_coeff=zeros(num_valid_sc,1);
for i=1:num_valid_sc
    if i==1
        chan_coeff(i)=(chan_coeff_matrix(1)+chan_coeff_matrix(2))/2;
    elseif i==52
        chan_coeff(i)=(chan_coeff_matrix(51)+chan_coeff_matrix(52))/2;
    else
        chan_coeff(i)=(chan_coeff_matrix(i-1)+chan_coeff_matrix(i)+chan_coeff_matrix(i+1))/3;
    end
end

chan_coeff=repmat(chan_coeff,1,20);
data_freq2=data_freq.*chan_coeff;
data_f_frm2 = zeros(num_sc, num_symbol_per_frame);
data_f_frm2(POS_VALID_SC,:) = data_freq2./phase_corr;
f_tx2=data_f_frm2;
intf=ifft(f_tx2)*sqrt(num_sc);
intf=[intf(end-15:end,:); intf];

load time_signal_tx1

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Interference Subtraction
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subtraction from original received signal
signal=ref_signal-alph*intf;
rx_signal_frame=reshape(signal,1,[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Process on the remaining desired signal %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % received frequency frame
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_signal_one_frame = rx_signal_frame;
freq_signal_one_frame = zeros(num_ant_ap, num_valid_sc, num_symbol_per_frame);
for rx_idx = 1:num_ant_ap
    time_signal_tmp1 = time_signal_one_frame(rx_idx,:);
    time_signal_tmp2 = reshape(time_signal_tmp1, num_sc+16, num_symbol_per_frame);
    time_signal_tmp3 = time_signal_tmp2(17:end,:);
    freq_signal_tmp1 = fft(time_signal_tmp3)/sqrt(num_sc);
    freq_signal_tmp2 = freq_signal_tmp1(POS_VALID_SC,:);
    freq_signal_one_frame(rx_idx,:,:) = freq_signal_tmp2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute EYY and EYX
EYY_SUM = zeros(num_ant_ap, num_ant_ap, num_valid_sc);
EYX_SUM = zeros(num_ant_ap, num_valid_sc);
ltf_dat=(1-alph)*LTF_DAT_GRP(1,:);
stf_dat=(1-alph)*STF_DAT_GRP(1,:);
for sc_idx = 1:num_valid_sc
    % TSF (OFDM 1)
        Y = freq_signal_one_frame(:, sc_idx, 1);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = stf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % TSF (OFDM 2)
        Y = freq_signal_one_frame(:, sc_idx, 2);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = stf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % LTF (OFDM 3)
        Y = freq_signal_one_frame(:, sc_idx, 3);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = ltf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
        % LTF (OFDM 4)
        Y = freq_signal_one_frame(:, sc_idx, 4);
        EYY_SUM(:,:,sc_idx) = squeeze(EYY_SUM(:,:,sc_idx)) + Y*Y';
        X = ltf_dat(1, sc_idx);
        EYX_SUM(:,sc_idx) = squeeze(EYX_SUM(:,sc_idx)) + Y*X';
end
    % compute G
    G_ARR_PER_FRAME = zeros(num_ant_ap, num_valid_sc);
for sc_idx = 1:num_valid_sc
        EYY = zeros(num_ant_ap, num_ant_ap);
        EYX = zeros(num_ant_ap, 1);
        sc_lwbd = max(1, sc_idx-num_averaging_sc);
        sc_upbd = min(sc_idx+num_averaging_sc, num_valid_sc);
        for kk = sc_lwbd:sc_upbd
            EYY = EYY + squeeze(EYY_SUM(:,:,kk));
            EYX = EYX + squeeze(EYX_SUM(:,kk));
        end
        G = pinv(EYY)*(EYX);
        G_ARR_PER_FRAME(:,sc_idx) = G;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decoded data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
decoded_data_frame = zeros(num_valid_sc,num_symbol_per_frame);
for sc_idx = 1:num_valid_sc
        G = G_ARR_PER_FRAME(:, sc_idx);
        Y = squeeze(freq_signal_one_frame(:, sc_idx, :));
        decoded_data_frame(sc_idx,:) = G'*Y;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decoded data with phase compension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bypass_freq_synch==0
PILOT_DAT_GRP=PILOT_DAT_GRP;   
phase_offset1 = decoded_data_frame(POS_PILOT_SC(1),5:end)./PILOT_DAT_GRP(1,5:end);
phase_offset2 = decoded_data_frame(POS_PILOT_SC(2),5:end)./PILOT_DAT_GRP(2,5:end);
phase_offset3 = decoded_data_frame(POS_PILOT_SC(3),5:end)./PILOT_DAT_GRP(3,5:end);
phase_offset4 = decoded_data_frame(POS_PILOT_SC(4),5:end)./PILOT_DAT_GRP(4,5:end);

phase_offset=ones(52,20);
phase_offset(1:13,5:end) =repmat(phase_offset1,13,1);
phase_offset(14:27,5:end)=repmat(phase_offset2,14,1);
phase_offset(28:40,5:end)=repmat(phase_offset3,13,1);
phase_offset(41:52,5:end)=repmat(phase_offset4,12,1);
decoded_data_frame=decoded_data_frame.*phase_corr;
else
    disp('1st Phase correction is bypassed')
end

decoded_payload_w_phase_comp = decoded_data_frame(POS_PAYLOAD_SC,5:num_symbol_per_frame);  
decoded_payload_w_phase_comp_arr = decoded_payload_w_phase_comp(:);
rx_data=decoded_payload_w_phase_comp_arr;
rx_data=rx_data/(1-alph);
demoded_data=qamdemod(rx_data,M,'UnitAveragePower',true); 
ref_points=qammod(demoded_data,M,'UnitAveragePower',true); 
err=abs(rx_data-ref_points);
evm=10*log10(mean(err.^2));
evm_per_subc=zeros(48,1);
ref1=reshape(ref_points,48,[]);
data1=reshape(rx_data,48,[]);
err=abs(ref1-data1);
evm_per_subc=10*log10(mean(err'.^2));
avm=sort(evm_per_subc);
fprintf('Our approach: %f',evm)
fprintf('Avg EVM on all: %f, Evm over top 24 subc: %f \n Best subc: %f, Worst subc: %f\n\n',evm, sum(avm(1:24))/24,avm(1),avm(end))
h = figure(); 
hold on; grid on;box on;
title(['Constellation of desired signal at strong user EVM=',num2str(evm),'dB']);            
scatter(real(decoded_payload_w_phase_comp_arr),imag(decoded_payload_w_phase_comp_arr),'ko','linewidth',0.5,'MarkerFaceColor','g');
axis([-2, 2, -2, 2])


h = figure(); 
hold on; grid on;box on;
title(['Re-scaled constellation of desired signal EVM=',num2str(evm)]);            
scatter(real(decoded_payload_w_phase_comp_arr)/(1-alph),imag(decoded_payload_w_phase_comp_arr)/(1-alph),'ko','linewidth',0.5,'MarkerFaceColor','g');
axis([-2, 2, -2, 2])