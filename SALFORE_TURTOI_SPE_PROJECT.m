clear all; clc;
close all;

%% =====================Input Parameters Definition======================%%
N = 256;                                    % The size of the transmitter IFFT and the receiver FFT
M = 8;                                     % Symbols per block
Q = N/M;                                    % The bandwidth spreading factor for IFDMA
Q_D = Q-1;                                  % The bandwidth spreeding factor for a DFDMA
cp = 32;                                    % CP length.
SP.subband = 2;                             % Set the subband location
order = 16;                                 % Modulation orders
modulation_types = 'qam';                   % Modulation types: QAM or PSK
equalizer_type = 'MMSE';                    % Equalization type: MMSE or ZFEQ
SNR = 0:1:30;                               % Simulated SNR range
N_simu = 5;                                % The number of simulations
alpha = 0.2;                                % Non-orthogonality parameter
channel = [0.96 10^(-9.7/20) 10^(-22.8/20)];   % Channels based on 3GPP TS 25.104.
channel = channel/sqrt(sum(channel.^2));    % Normalize the channel.
num_users = N / M;                          % Number of users in the system
H = fft(channel,N);                         % Frequency expression of the channel response

%% Transmitter
%=================== symbol error initialization =========================%
Ifdma_SER = zeros(1,length(SNR));
Dfdma_SER = zeros(1,length(SNR));
Lfdma_SER = zeros(1,length(SNR));
Nfdma_SER = zeros(1,length(SNR));
paprI = zeros (1,N_simu);
paprD = zeros (1,N_simu);
paprL = zeros (1,N_simu);
paprN = zeros (1,N_simu);
figure % Constellation figure
hold on;
for n = 1:length(SNR)
    % Initialize the error count.
    Ifdma_error_number = zeros(1, num_users);
    Dfdma_error_number = zeros(1, num_users);
    Lfdma_error_number = zeros(1, num_users);
    Nfdma_error_number = zeros(1, num_users);
    for u = 1:num_users
        for k = 1:N_simu
            %=================== Generate random data block===========================%
            input_symbol = randi([0, (order-1)], 1, M);
            %==============================Modulation===============================%
            if strcmp(modulation_types, 'qam')
                input_signal = qammod(input_symbol, order);
            elseif strcmp(modulation_types, 'psk')
                input_signal = pskmod(input_symbol, order);
            end
            %========================Tranform to frequency domain=====================%
            input_signal_fft = fft(input_signal, M); 
            %===========================Subcarrier Mapping============================%
            % Initialize subcarrier mappings
            Ifdma_mapping = zeros(1, N);
            Dfdma_mapping = zeros(1, N);
            Lfdma_mapping = zeros(1, N);
            Nfdma_mapping = zeros(1, N);

            % Applying Mapping
            Ifdma_mapping(1+SP.subband:Q:N) = input_signal_fft;
            if Q == 1
                Dfdma_mapping = Ifdma_mapping;
            else
                Dfdma_mapping(1+SP.subband:Q_D:Q_D*M) = input_signal_fft;
            end
            Lfdma_mapping(1:M+M*SP.subband) = input_signal_fft;
            
            nonOrthogonalMapping = randperm(N, M); % Randomly select M subcarriers
            Nfdma_mapping(nonOrthogonalMapping) = (1 - alpha) * input_signal_fft;

            %========================Tranform back to time domain=====================%
            Ifdma_IFFT = ifft(Ifdma_mapping, N);
            Dfdma_IFFT = ifft(Dfdma_mapping, N);
            Lfdma_IFFT = ifft(Lfdma_mapping, N);
            Nfdma_IFFT = ifft(Nfdma_mapping, N);
            %==========================Add a cyclic prefix============================%
            % Cyclic Prefix Addition
            Ifdma_cyclic = [Ifdma_IFFT(N-cp+1:N) Ifdma_IFFT];
            Dfdma_cyclic = [Dfdma_IFFT(N-cp+1:N) Dfdma_IFFT];
            Lfdma_cyclic = [Lfdma_IFFT(N-cp+1:N) Lfdma_IFFT];
            Nfdma_cyclic = [Nfdma_IFFT(N-cp+1:N) Nfdma_IFFT];

            %% CHANNEL
            % Multi-path channel
            Ifdma_signal = filter(channel, 1, Ifdma_cyclic);
            Dfdma_signal = filter(channel, 1, Dfdma_cyclic);
            Lfdma_signal = filter(channel, 1, Lfdma_cyclic);
            Nfdma_signal = filter(channel, 1, Nfdma_cyclic);

            % Generate AWGN
            Noise  = (randn(1, N+cp) + 1i*randn(1, N+cp))/sqrt(2); % N+cp by considering the added cyclic symbols
            noisePower = 10^(-SNR(n)/10);

            % Add AWGN to the transmitted signal
            Ifdma_signal = Ifdma_signal + sqrt(noisePower/M) * Noise;
            Dfdma_signal = Dfdma_signal + sqrt(noisePower/M) * Noise;
            Lfdma_signal = Lfdma_signal + sqrt(noisePower/M) * Noise;
            Nfdma_signal = Nfdma_signal + sqrt(noisePower/M) * Noise;

            %% Receiver
            %=========================Remove the cyclic prefix========================%
            % Removing cyclic prefix
            Ifdma_received = Ifdma_signal(cp+1:N+cp);
            Dfdma_received = Dfdma_signal(cp+1:N+cp);
            Lfdma_received = Lfdma_signal(cp+1:N+cp);
            Nfdma_received = Nfdma_signal(cp+1:N+cp);

            %========================Tranform to frequency domain=====================%
            % Applying FFT Operation
            Ifdma_received = fft(Ifdma_received, N);
            Dfdma_received = fft(Dfdma_received, N);
            Lfdma_received = fft(Lfdma_received, N);
            Nfdma_received = fft(Nfdma_received, N);

            %===========================Subcarrier De-Mapping=========================%
            % Applying demapping
            Ifdma_received = Ifdma_received(1+SP.subband:Q:N);
            if Q == 1
                Dfdma_received = Ifdma_received;
            else
                Dfdma_received = Dfdma_received(1+SP.subband:Q_D:Q_D*M);
            end
            Lfdma_received = Lfdma_received(1:M+M*SP.subband);

            Nfdma_received = Nfdma_received(nonOrthogonalMapping);

            %====================Perform frequency equalization=======================%
            % Channel response of the subcarriers
            H_Ifdma = H(1+SP.subband:Q:N);
            if Q == 1
                H_Dfdma =  H_Ifdma;
            else
                H_Dfdma = H(1+SP.subband:Q_D:Q_D*M);
            end
            H_Lfdma = H(1:M+M*SP.subband);
            H_Nfdma = H(nonOrthogonalMapping);

            if strcmp(equalizer_type, 'ZFEQ')
                Ifdma_received = Ifdma_received ./ H_Ifdma;
                Dfdma_received = Dfdma_received ./ H_Dfdma;
                Lfdma_received = Lfdma_received ./ H_Lfdma;
                Nfdma_received = Nfdma_received ./ H_Nfdma;

            elseif strcmp(equalizer_type, 'MMSE')
                C_Ifdma = conj(H_Ifdma)./(conj(H_Ifdma).*H_Ifdma + 10^(-SNR(n)/10));
                C_Dfdma = conj(H_Dfdma)./(conj(H_Dfdma).*H_Dfdma + 10^(-SNR(n)/10));
                C_Lfdma = conj(H_Lfdma)./(conj(H_Lfdma).*H_Lfdma + 10^(-SNR(n)/10));
                C_Nfdma = conj(H_Nfdma)./(conj(H_Nfdma).*H_Nfdma + 10^(-SNR(n)/10));
                Ifdma_received  = Ifdma_received .* C_Ifdma;
                Dfdma_received  = Dfdma_received .* C_Dfdma;
                Lfdma_received  = Lfdma_received .* C_Lfdma;
                Nfdma_received  = Nfdma_received .* C_Nfdma;
            end
            %========================Tranform back to time domain=====================%
            Ifdma_received = ifft(Ifdma_received, M);
            Dfdma_received = ifft(Dfdma_received, M);
            Lfdma_received = ifft(Lfdma_received, M);
            Nfdma_received = ifft(Nfdma_received, M);
            %=============================De-Modulation===============================%
            % Symbol detection
            if strcmp(modulation_types, 'qam')
                Ifdma_symbol = qamdemod(Ifdma_received, order);
                Dfdma_symbol = qamdemod(Dfdma_received, order);
                Lfdma_symbol = qamdemod(Lfdma_received, order);
                Nfdma_symbol = qamdemod(Nfdma_received, order);
            elseif strcmp(modulation_types, 'psk')
                Ifdma_symbol = pskdemod(Ifdma_received, order);
                Dfdma_symbol = pskdemod(Dfdma_received, order);
                Lfdma_symbol = pskdemod(Lfdma_received, order);
                Nfdma_symbol = pskdemod(Nfdma_received, order);
            end

            %% PAPR Calculation
            paprI(k) = 10*log10(max(abs(Ifdma_received).^2)/mean(abs(Ifdma_received).^2));
            paprD(k) = 10*log10(max(abs(Dfdma_received).^2)/mean(abs(Dfdma_received).^2));
            paprL(k) = 10*log10(max(abs(Lfdma_received).^2)/mean(abs(Lfdma_received).^2));
            paprN(k) = 10*log10(max(abs(Nfdma_received).^2)/mean(abs(Nfdma_received).^2));
            %% Error Calculation
            % Number of correctly received symbols
            Ifdma_correct_symbol = find((input_symbol - Ifdma_symbol) == 0);
            Dfdma_correct_symbol = find((input_symbol - Dfdma_symbol) == 0);
            Lfdma_correct_symbol = find((input_symbol - Lfdma_symbol) == 0);
            Nfdma_correct_symbol = find((input_symbol - Nfdma_symbol) == 0);
            % The number of errors
            Ifdma_error_number(u) = Ifdma_error_number(u) + (M - length(Ifdma_correct_symbol));
            Dfdma_error_number(u) = Dfdma_error_number(u) + (M - length(Dfdma_correct_symbol));
            Lfdma_error_number(u) = Lfdma_error_number(u) + (M - length(Lfdma_correct_symbol));
            Nfdma_error_number(u) = Nfdma_error_number(u) + (M - length(Nfdma_correct_symbol));
        end
        % Calculate the SER
        Ifdma_SER(n) = Ifdma_error_number(u) / (M * N_simu);
        Dfdma_SER(n) = Dfdma_error_number(u) / (M * N_simu);
        Lfdma_SER(n) = Lfdma_error_number(u) / (M * N_simu);
        Nfdma_SER(n) = Nfdma_error_number(u) / (M * N_simu);
    end
    % Plot the received symbols on the constellation diagram
    plot(real(Ifdma_received), imag(Ifdma_received), 'bo');
    plot(real(Dfdma_received), imag(Dfdma_received), 'bo');
    plot(real(Lfdma_received), imag(Lfdma_received), 'bo');
    plot(real(Nfdma_received), imag(Nfdma_received), 'bo');
end
title('Constellation Diagram');
xlabel('In-Phase');
ylabel('Quadrature');
grid on;
hold off;

%% Plotting
% Plot SER curves
figure;
% Find first zero symbol error rate
first_zero1 = find(Ifdma_SER(:) == 0, 1, 'first');
first_zero2 = find(Dfdma_SER(:) == 0, 1, 'first');
first_zero3 = find(Lfdma_SER(:) == 0, 1, 'first');
first_zero4 = find(Nfdma_SER(:) == 0, 1, 'first');
% Plot symbol error rate
semilogy(SNR(1:first_zero1), max(10^-8, Ifdma_SER(1:first_zero1)), 'LineWidth', 2);
hold on;
semilogy(SNR(1:first_zero2), max(10^-8, Dfdma_SER(1:first_zero2)), 'LineWidth', 2);
semilogy(SNR(1:first_zero3), max(10^-8, Lfdma_SER(1:first_zero3)), 'LineWidth', 2);
semilogy(SNR(1:first_zero4), max(10^-8, Nfdma_SER(1:first_zero4)), 'LineWidth', 2);
hold off;
legend('IFDMA', 'DFDMA', 'LFDMA', 'NFDMA', 'FontSize', 10);
title(['SER - ','Block Length ', num2str(M),' ',upper(modulation_types), ' Order = ', num2str(order), ' Eq Type = ', equalizer_type]);
xlabel('Signal to Noise Ratio in [dB]', 'FontSize', 10);
ylabel('Symbol Error rate', 'FontSize', 10);
set(gca, 'FontSize', 12);
axis([0 30 10^-8 1]);

% Calculate the ECDF for PAPR values
sorted_paprI = sort(paprI);
sorted_paprD = sort(paprD);
sorted_paprL = sort(paprL);
sorted_paprN = sort(paprN);
ecdf = (1:N_simu) / N_simu;
min_papr = min([min(sorted_paprI), min(sorted_paprD), min(sorted_paprL), min(sorted_paprN)]);
max_papr = max([max(sorted_paprI), max(sorted_paprD), max(sorted_paprL), max(sorted_paprN)]);

% Plot the ECDF and quartiles
figure;
plot(sorted_paprI, ecdf, 'LineWidth', 2, 'DisplayName', 'IFDMA');
hold on;
plot(sorted_paprD, ecdf, 'LineWidth', 2, 'DisplayName', 'DFDMA');
plot(sorted_paprL, ecdf, 'LineWidth', 2, 'DisplayName', 'LFDMA');
plot(sorted_paprN, ecdf, 'LineWidth', 2, 'DisplayName', 'NFDMA');

% Calculate and plot the first and third quartiles
first_quartile_I = quantile(paprI, 0.25);
third_quartile_I = quantile(paprI, 0.75);
first_quartile_D = quantile(paprD, 0.25);
third_quartile_D = quantile(paprD, 0.75);
first_quartile_L = quantile(paprL, 0.25);
third_quartile_L = quantile(paprL, 0.75);
first_quartile_N = quantile(paprN, 0.25);
third_quartile_N = quantile(paprN, 0.75);

plot([first_quartile_I, first_quartile_I], [0, 0.25], 'r--', 'LineWidth', 0.5);
plot([third_quartile_I, third_quartile_I], [0, 0.75], 'r--', 'LineWidth', 0.5);
plot([first_quartile_D, first_quartile_D], [0, 0.25], 'g--', 'LineWidth', 0.5);
plot([third_quartile_D, third_quartile_D], [0, 0.75], 'g--', 'LineWidth', 0.5);
plot([first_quartile_L, first_quartile_L], [0, 0.25], 'b--', 'LineWidth', 0.5);
plot([third_quartile_L, third_quartile_L], [0, 0.75], 'b--', 'LineWidth', 0.5);
plot([first_quartile_N, first_quartile_N], [0, 0.25], 'm--', 'LineWidth', 0.5);
plot([third_quartile_N, third_quartile_N], [0, 0.75], 'm--', 'LineWidth', 0.5);

plot([min_papr, max_papr], [0.25, 0.25], 'k--', 'LineWidth', 1.5);
plot([min_papr, max_papr], [0.75, 0.75], 'k--', 'LineWidth', 1.5);

hold off;
title(['PAPR '  ' '  ' ' equalizer_type], 'FontSize', 12);
xlabel('PAPR (dB)', 'FontSize', 12);
ylabel('Cumulative Probability', 'FontSize', 12);
legend('IFDMA', 'DFDMA', 'LFDMA', 'NFDMA', 'Location', 'southeast', 'FontSize', 10);
xlim([min_papr, max_papr]);

% Create box plot for PAPR values
figure;
boxplot([paprI', paprD', paprL', paprN'], 'Labels', {'IFDMA', 'DFDMA', 'LFDMA', 'NFDMA'});
title('Box Plot of PAPR Values', 'FontSize', 12);
xlabel('Subcarrier Mapping Scheme', 'FontSize', 10);
ylabel('PAPR (dB)', 'FontSize', 10);
