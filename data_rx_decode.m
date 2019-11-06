clear all;
clc;
close all;

qpsk_symbols=[1/sqrt(2)+1j/sqrt(2) 1/sqrt(2)-1j/sqrt(2) -1/sqrt(2)+1j/sqrt(2) -1/sqrt(2)-1j/sqrt(2)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_debug_enabled = 0;
num_ant_ap = 1;
num_stm_ap = 1;
num_sample_processing = 80000;
num_sc = 64;
num_symbol_per_frame = 100;
num_sample_shift = 1e6;
num_valid_sc = 52;
const_correl_threhold = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STF 
STF_POS = [3  7  11  15  19  23  30  34  38  42  46  50];
STF_DAT_GRP = zeros(1, 52);    
STF_DAT_GRP(1,STF_POS) = sqrt(2)*[-1-1i  -1+1i  -1-1i  -1-1i  -1-1i  1-1i  -1-1i  -1-1i  1+1i  1-1i  1+1i  -1+1i];

%%% LTF
LTF_DAT_GRP = zeros(3, 52);
LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(2,:) = [1 1 -1 -1 1 1 -1 -1 1 1 -1 -1 -1 1 -1 1 1 1 -1 -1 -1 -1 1 1 1 -1 1 1 1 -1 -1 -1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 -1 -1 1 -1 -1 1 1 -1];
LTF_DAT_GRP(3,:) = [1 -1 -1 -1 -1 1 -1 -1 -1 1 1 -1 1 1 1 1 1 -1 -1 1 -1 -1 -1 -1 -1 1 -1 1 -1 -1 1 1 -1 -1 -1 1 1 -1 1 -1 1 -1 1 1 -1 -1 1 1 1 -1 1 -1];

%%%% pilot
POS_VALID_SC = [39:64 2:27];
POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
POS_PILOT_SC = [6 20 33 47];
PILOT_DAT_GRP=ones(4,num_symbol_per_frame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data from user 1
file_index=1;
file_id = fopen(['C:\Users\pedra\Desktop\NOMA stuff_2nd\NOMA_2X3\signals\noma_rx_u' num2str(file_index) '.dat'], 'rb');%
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
ltf_freq_signal(POS_VALID_SC) = LTF_DAT_GRP(1,:);
lte_time_signal = ifft(ltf_freq_signal)*sqrt(num_sc);
ltf_waveform = [lte_time_signal(end-15:end) lte_time_signal lte_time_signal(end-15:end) lte_time_signal];
cross_correl_value = zeros(1,num_sample_processing);
for ii = 161:num_sample_processing-160
    signal_segment = rx_signal_buf(ii:ii+160-1);
    cross_correl_value(ii) =  ltf_waveform*signal_segment'/(norm(ltf_waveform)*norm(signal_segment));
end
[max_correl_value, max_correl_pos] = max(abs(cross_correl_value(1:end-num_symbol_per_frame*80)));
begofframe = max_correl_pos;
if (is_debug_enabled == 11 && max_correl_value > const_correl_threhold)
    figure(); hold on; grid on; box on;
    plot(abs(cross_correl_value),'k-');
    xlabel('sample index');
    ylabel('correlation');
    axis([0 length(rx_signal_buf) 0 1]);
    title('Cross correlation of received signal at user 1');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the begining of a frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_signal_buf=rx_signal_buf(begofframe:end);
rx_signal_frame=rx_signal_buf(1:num_symbol_per_frame*80);
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
if is_debug_enabled == 11
   disp(['cp-based freq offset estimate: ' num2str(delta_phase_per_sample_cp_based)]);
   disp(['ltf-based freq offset estimate: ' num2str(delta_phase_per_sample_ltf_based)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% received frequency frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_signal_one_frame = rx_signal_frame;
time_signal_tmp1 = time_signal_one_frame;
time_signal_tmp2 = reshape(time_signal_tmp1, num_sc+16, num_symbol_per_frame);
time_signal_tmp3 = time_signal_tmp2(17:end,:);
freq_signal_tmp1 = fft(time_signal_tmp3)/sqrt(num_sc);
freq_signal_tmp2 = freq_signal_tmp1(POS_VALID_SC,:);
freq_signal_one_frame= freq_signal_tmp2;
if is_debug_enabled==11
   figure()
   freq_signal_one_frame=freq_signal_one_frame(1:4,:);
   x=freq_signal_one_frame(:).';
   scatter(real(x),imag(x),'ko');
   title('Received signal at user 1-before channel equalization')
end
%% zf
h=0.5*(freq_signal_tmp2(:,1)+freq_signal_tmp2(:,2))./(LTF_DAT_GRP(1,:)).';
h2=repmat(h,1,num_symbol_per_frame);
freq_signal_tmp2=freq_signal_tmp2./h2;
freq_signal_tmp3=freq_signal_tmp2(POS_PAYLOAD_SC(1:4),7:end);
freq_signal_tmp3=freq_signal_tmp3(:).';
if is_debug_enabled==11
   figure()
   scatter(real(freq_signal_tmp3),imag(freq_signal_tmp3),'ko');
   hold on
   scatter(real(qpsk_symbols),imag(qpsk_symbols),'ro','linewidth',2);
   legend('Received data','QPSK symbols')
   title('Equlized signal at user 1-before phase correction')
end
rec_pilots=freq_signal_tmp2(POS_PILOT_SC,7:end);
phi=angle(rec_pilots(1,:));
delta_phi=zeros(4,length(phi(1,:))-1);
for i=2:length(phi(1,:))
    delta_phi(:,i-1)=phi(:,i)-phi(:,i-1);
end
if is_debug_enabled==11
   figure()
   plot(1:length(phi(1,:))-1,delta_phi(1,:),'-ob')
   hold on
   plot(1:length(phi(1,:))-1,delta_phi(2,:),'-sr')
   plot(1:length(phi(1,:))-1,delta_phi(1,:),'-vg')
   plot(1:length(phi(1,:))-1,delta_phi(1,:),'-*m')
   legend('Pilot 1','Pilot 2','Pilot 3','Pilot 4')
   title('\Delta\phi at user 1')
end

phase_offset=zeros(4,1);
for i=1:4;
    phase_offset(i)=mean(delta_phi(i,:));
end
decoded_data_frame=zeros(52,num_symbol_per_frame-6);
decoded_data_frame(1:13,:) = exp(-1j*phase_offset(1))*freq_signal_tmp2(1:13,7:end);
decoded_data_frame(14:27,:) = exp(-1j*phase_offset(2))*freq_signal_tmp2(14:27,7:end);
decoded_data_frame(28:40,:) = exp(-1j*phase_offset(3))*freq_signal_tmp2(28:40,7:end);
decoded_data_frame(41:52,:) = exp(-1j*phase_offset(4))*freq_signal_tmp2(41:52,7:end);
x=decoded_data_frame(1:4,:);
x=x(:).';
figure()
scatter(real(x),imag(x),'ko');
hold on
scatter(real(qpsk_symbols),imag(qpsk_symbols),'ro','linewidth',2);
legend('Received data','QPSK symbols')
title('Desired signal at the weakest user')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data from user 2
file_index=2;
file_id = fopen(['C:\Users\pedra\Desktop\NOMA_2X3\signals\noma_rx_u' num2str(file_index) '.dat'], 'rb');%
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
ltf_freq_signal(POS_VALID_SC) = LTF_DAT_GRP(1,:);
lte_time_signal = ifft(ltf_freq_signal)*sqrt(num_sc);
ltf_waveform = [lte_time_signal(end-15:end) lte_time_signal lte_time_signal(end-15:end) lte_time_signal];
cross_correl_value = zeros(1,num_sample_processing);
for ii = 161:num_sample_processing-160
    signal_segment = rx_signal_buf(ii:ii+160-1);
    cross_correl_value(ii) =  ltf_waveform*signal_segment'/(norm(ltf_waveform)*norm(signal_segment));
end
[max_correl_value, max_correl_pos] = max(abs(cross_correl_value(1:end-num_symbol_per_frame*80)));
begofframe = max_correl_pos;
if (is_debug_enabled == 11 && max_correl_value > const_correl_threhold)
    figure(); hold on; grid on; box on;
    plot(abs(cross_correl_value),'k-');
    xlabel('sample index');
    ylabel('correlation');
    axis([0 length(rx_signal_buf) 0 1]);
    title('Cross correlation of received signal at user 2');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the begining of a frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_signal_buf=rx_signal_buf(begofframe:end);
rx_signal_frame=rx_signal_buf(1:num_symbol_per_frame*80);
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
if is_debug_enabled == 11
   disp(['cp-based freq offset estimate: ' num2str(delta_phase_per_sample_cp_based)]);
   disp(['ltf-based freq offset estimate: ' num2str(delta_phase_per_sample_ltf_based)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% received frequency frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_signal_one_frame = rx_signal_frame;
time_signal_tmp1 = time_signal_one_frame;
time_signal_tmp2 = reshape(time_signal_tmp1, num_sc+16, num_symbol_per_frame);
time_signal_tmp3 = time_signal_tmp2(17:end,:);
freq_signal_tmp1 = fft(time_signal_tmp3)/sqrt(num_sc);
freq_signal_tmp2 = freq_signal_tmp1(POS_VALID_SC,:);
freq_signal_one_frame= freq_signal_tmp2;
if is_debug_enabled==11
   figure()
   freq_signal_one_frame=freq_signal_one_frame(1:4,:);
   x=freq_signal_one_frame(:).';
   scatter(real(x),imag(x),'ko');
   title('Received signal at user 2-before channel equalization')
end
%% zf
ref_sig=freq_signal_tmp2;
h=0.5*(freq_signal_tmp2(:,1)+freq_signal_tmp2(:,2))./(LTF_DAT_GRP(1,:)).';
h2=repmat(h,1,num_symbol_per_frame-6);
freq_signal_tmp2(:,7:end)=freq_signal_tmp2(:,7:end)./h2;
freq_signal_tmp3=freq_signal_tmp2(POS_PAYLOAD_SC(1:4),7:end);
freq_signal_tmp3=freq_signal_tmp3(:).';
if is_debug_enabled==11
   figure()
   scatter(real(freq_signal_tmp3),imag(freq_signal_tmp3),'ko');
   hold on
   scatter(real(qpsk_symbols),imag(qpsk_symbols),'ro','linewidth',2);
   legend('Received data','QPSK symbols')
   title('Equlized signal at user 2-before phase correction')
end
rec_pilots=freq_signal_tmp2(POS_PILOT_SC,7:end);
phi=angle(rec_pilots(1,:));
delta_phi=zeros(4,length(phi(1,:))-1);
for i=2:length(phi(1,:))
    delta_phi(:,i-1)=phi(:,i)-phi(:,i-1);
end
if is_debug_enabled==11
   figure()
   plot(1:length(phi(1,:))-1,delta_phi(1,:),'-ob')
   hold on
   plot(1:length(phi(1,:))-1,delta_phi(2,:),'-sr')
   plot(1:length(phi(1,:))-1,delta_phi(1,:),'-vg')
   plot(1:length(phi(1,:))-1,delta_phi(1,:),'-*m')
   legend('Pilot 1','Pilot 2','Pilot 3','Pilot 4')
   title('\Delta\phi at user 1')
end

phase_offset=zeros(4,1);
for i=1:4;
    phase_offset(i)=mean(delta_phi(i,:));
end
decoded_data_frame=zeros(52,num_symbol_per_frame);
decoded_data_frame(1:13,:) = exp(-1j*phase_offset(1))*freq_signal_tmp2(1:13,1:end);
decoded_data_frame(14:27,:) = exp(-1j*phase_offset(2))*freq_signal_tmp2(14:27,1:end);
decoded_data_frame(28:40,:) = exp(-1j*phase_offset(3))*freq_signal_tmp2(28:40,1:end);
decoded_data_frame(41:52,:) = exp(-1j*phase_offset(4))*freq_signal_tmp2(41:52,1:end);

%%%%%%%%% I'm not sure %%%%%%%%%%%%%%%
ref_sig(1:13,:) = exp(-1j*phase_offset(1))*ref_sig(1:13,1:end);
ref_sig(14:27,:) = exp(-1j*phase_offset(2))*ref_sig(14:27,1:end);
ref_sig(28:40,:) = exp(-1j*phase_offset(3))*ref_sig(28:40,1:end);
ref_sig(41:52,:) = exp(-1j*phase_offset(4))*ref_sig(41:52,1:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=decoded_data_frame(1:4,7:end);
x2=x(:).';
figure()
scatter(real(x2),imag(x2),'ko');
hold on
scatter(real(qpsk_symbols),imag(qpsk_symbols),'ro','linewidth',2);
legend('Received data','QPSK symbols')
title('Received signal at user 2')
% SIC
x_hat=qammod(qamdemod(x,4,'UnitAveragePower',true),4,'UnitAveragePower',true);
% interference frame reconstrunction
Intf=zeros(52,num_symbol_per_frame);
Intf(1:4,7:end)=x_hat;
ref_sig1=ref_sig-repmat(h,1,num_symbol_per_frame).*Intf;
h=0.5*(ref_sig1(:,3)+ref_sig1(:,4))./(LTF_DAT_GRP(2,:)).';
h2=repmat(h,1,num_symbol_per_frame-6);
ref_sig1(:,7:end)=ref_sig1(:,7:end)./h2;
x=ref_sig1(1:4,7:end);
x=x(:).';
figure()
scatter(real(x),imag(x),'ko');
hold on
scatter(real(qpsk_symbols),imag(qpsk_symbols),'ro','linewidth',2);
legend('Received data','QPSK symbols')
title('Desired signal at user 2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data from user 3
file_index=3;
file_id = fopen(['C:\Users\pedra\Desktop\NOMA_2X3\signals\noma_rx_u' num2str(file_index) '.dat'], 'rb');%
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
ltf_freq_signal(POS_VALID_SC) = LTF_DAT_GRP(1,:);
lte_time_signal = ifft(ltf_freq_signal)*sqrt(num_sc);
ltf_waveform = [lte_time_signal(end-15:end) lte_time_signal lte_time_signal(end-15:end) lte_time_signal];
cross_correl_value = zeros(1,num_sample_processing);
for ii = 161:num_sample_processing-160
    signal_segment = rx_signal_buf(ii:ii+160-1);
    cross_correl_value(ii) =  ltf_waveform*signal_segment'/(norm(ltf_waveform)*norm(signal_segment));
end
[max_correl_value, max_correl_pos] = max(abs(cross_correl_value(1:end-num_symbol_per_frame*80)));
begofframe = max_correl_pos;
if (is_debug_enabled == 11 && max_correl_value > const_correl_threhold)
    figure(); hold on; grid on; box on;
    plot(abs(cross_correl_value),'k-');
    xlabel('sample index');
    ylabel('correlation');
    axis([0 length(rx_signal_buf) 0 1]);
    title('Cross correlation of received signal at user 2');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the begining of a frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_signal_buf=rx_signal_buf(begofframe:end);
rx_signal_frame=rx_signal_buf(1:num_symbol_per_frame*80);
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
if is_debug_enabled == 11
   disp(['cp-based freq offset estimate: ' num2str(delta_phase_per_sample_cp_based)]);
   disp(['ltf-based freq offset estimate: ' num2str(delta_phase_per_sample_ltf_based)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% received frequency frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_signal_one_frame = rx_signal_frame;
time_signal_tmp1 = time_signal_one_frame;
time_signal_tmp2 = reshape(time_signal_tmp1, num_sc+16, num_symbol_per_frame);
time_signal_tmp3 = time_signal_tmp2(17:end,:);
freq_signal_tmp1 = fft(time_signal_tmp3)/sqrt(num_sc);
freq_signal_tmp2 = freq_signal_tmp1(POS_VALID_SC,:);
freq_signal_one_frame= freq_signal_tmp2;
if is_debug_enabled==11
   figure()
   freq_signal_one_frame=freq_signal_one_frame(1:4,:);
   x=freq_signal_one_frame(:).';
   scatter(real(x),imag(x),'ko');
   title('Received signal at user 3-before channel equalization')
end
%% zf
ref_sig=freq_signal_tmp2;
h=0.5*(freq_signal_tmp2(:,1)+freq_signal_tmp2(:,2))./(LTF_DAT_GRP(1,:)).';
h2=repmat(h,1,num_symbol_per_frame-6);
freq_signal_tmp2(:,7:end)=freq_signal_tmp2(:,7:end)./h2;
freq_signal_tmp3=freq_signal_tmp2(POS_PAYLOAD_SC(1:4),7:end);
freq_signal_tmp3=freq_signal_tmp3(:).';
if is_debug_enabled==11
   figure()
   scatter(real(freq_signal_tmp3),imag(freq_signal_tmp3),'ko');
   hold on
   scatter(real(qpsk_symbols),imag(qpsk_symbols),'ro','linewidth',2);
   legend('Received data','QPSK symbols')
   title('Equlized signal at user 3-before phase correction')
end
rec_pilots=freq_signal_tmp2(POS_PILOT_SC,7:end);
phi=angle(rec_pilots(1,:));
delta_phi=zeros(4,length(phi(1,:))-1);
for i=2:length(phi(1,:))
    delta_phi(:,i-1)=phi(:,i)-phi(:,i-1);
end
if is_debug_enabled==11
   figure()
   plot(1:length(phi(1,:))-1,delta_phi(1,:),'-ob')
   hold on
   plot(1:length(phi(1,:))-1,delta_phi(2,:),'-sr')
   plot(1:length(phi(1,:))-1,delta_phi(1,:),'-vg')
   plot(1:length(phi(1,:))-1,delta_phi(1,:),'-*m')
   legend('Pilot 1','Pilot 2','Pilot 3','Pilot 4')
   title('\Delta\phi at user 3')
end

phase_offset=zeros(4,1);
for i=1:4;
    phase_offset(i)=mean(delta_phi(i,:));
end
decoded_data_frame=zeros(52,num_symbol_per_frame);
decoded_data_frame(1:13,:) = exp(-1j*phase_offset(1))*freq_signal_tmp2(1:13,1:end);
decoded_data_frame(14:27,:) = exp(-1j*phase_offset(2))*freq_signal_tmp2(14:27,1:end);
decoded_data_frame(28:40,:) = exp(-1j*phase_offset(3))*freq_signal_tmp2(28:40,1:end);
decoded_data_frame(41:52,:) = exp(-1j*phase_offset(4))*freq_signal_tmp2(41:52,1:end);

%%%%%%%%% I'm not sure %%%%%%%%%%%%%%%
ref_sig(1:13,:) = exp(-1j*phase_offset(1))*ref_sig(1:13,1:end);
ref_sig(14:27,:) = exp(-1j*phase_offset(2))*ref_sig(14:27,1:end);
ref_sig(28:40,:) = exp(-1j*phase_offset(3))*ref_sig(28:40,1:end);
ref_sig(41:52,:) = exp(-1j*phase_offset(4))*ref_sig(41:52,1:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=decoded_data_frame(1:4,7:end);
x2=x(:).';
figure()
scatter(real(x2),imag(x2),'ko');
hold on
scatter(real(qpsk_symbols),imag(qpsk_symbols),'ro','linewidth',2);
legend('Received data','QPSK symbols')
title('Received signal at user 3')

% First iteration of SIC
x_hat=qammod(qamdemod(x,4,'UnitAveragePower',true),4,'UnitAveragePower',true);
Intf=zeros(52,num_symbol_per_frame);
Intf(1:4,7:end)=x_hat;
ref_sig=ref_sig-repmat(h,1,num_symbol_per_frame).*Intf;
h=0.5*(ref_sig(:,3)+ref_sig(:,4))./(LTF_DAT_GRP(2,:)).';
h2=repmat(h,1,num_symbol_per_frame-6);
ref_sig1=ref_sig(:,7:end)./h2;
x=ref_sig1(1:4,:);
x2=x(:).';
figure()
scatter(real(x2),imag(x2),'ko');
hold on
scatter(real(qpsk_symbols),imag(qpsk_symbols),'ro','linewidth',2);
legend('Received data','QPSK symbols')
title('Removing 1st interference at user 3');

% second iteration of SIC
x_hat=qammod(qamdemod(x,4,'UnitAveragePower',true),4,'UnitAveragePower',true);
Intf=zeros(52,num_symbol_per_frame);
Intf(1:4,7:end)=x_hat;
ref_sig=ref_sig-repmat(h,1,num_symbol_per_frame).*Intf;
h=0.5*(ref_sig(:,5)+ref_sig(:,6))./(LTF_DAT_GRP(3,:)).';
h2=repmat(h,1,num_symbol_per_frame-6);
ref_sig1=ref_sig(:,7:end)./h2;
x=ref_sig1(1:4,:);
x2=x(:).';
figure()
scatter(real(x2),imag(x2),'ko');
hold on
%scatter(real(qpsk_symbols),imag(qpsk_symbols),'ro','linewidth',2);
%legend('Received data','QPSK symbols')
title('Desired signal at user 3');

