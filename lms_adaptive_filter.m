% Variable-step-size LMS for adaptive noise cancellation.
%
% Compares standard LMS against two variants in which the step size mu(n)
% is driven by the instantaneous error, so the filter adapts quickly while
% the error is large and settles down as it shrinks.
%
% The desired signal is a sinusoid in additive white Gaussian noise; the
% filter input is that same signal delayed by one sample, making this a
% one-step linear predictor. Learning curves are averaged over 100
% independent runs at each of three SNRs.

close all
clear
clc

% --- Simulation parameters -------------------------------------------------
L        = 100;    % filter order (number of taps)
N        = 10000;  % samples per run
mu0      = 1e-5;   % fixed step size for standard LMS
alpha    = mu0;    % scaling for the variable step sizes
delta    = 0.01;   % error scale in the mu(n) laws
num_runs = 100;    % Monte Carlo runs to average over
lag      = 1;      % input delay, in samples

nn = 0:N-1;
x  = sin(2*pi*nn/20);          % input sinusoid, period 20 samples

SNR     = [10, 20, 30];        % dB
num_SNR = length(SNR);

% Mean-square error curves, accumulated across runs
e0_sr = zeros(N, num_SNR);     % standard LMS
e1_sr = zeros(N, num_SNR);     % variant 1
e2_sr = zeros(N, num_SNR);     % variant 2

for i = 1:num_SNR
    d(:, i)    = awgn(x, SNR(i), 'measured');           % noisy desired signal
    dlag(:, i) = [zeros(lag,1); d(1:end-lag, i)];       % delayed filter input

    for run = 1:num_runs
        W0 = zeros(L, N);
        W1 = zeros(L, N);
        W2 = zeros(L, N);

        mu1 = zeros(N, 1);
        mu2 = zeros(N, 1);

        y0 = zeros(N, 1);  e0 = zeros(N, 1);
        y1 = zeros(N, 1);  e1 = zeros(N, 1);
        y2 = zeros(N, 1);  e2 = zeros(N, 1);

        % --- Standard LMS: fixed step size ---------------------------------
        for n = L:N-1
            Xn     = dlag(n:-1:n-L+1, i);
            y0(n)  = Xn' * W0(:, n);
            e0(n)  = d(n, i) - y0(n);
            W0(:, n+1) = W0(:, n) + 2 * mu0 * e0(n) * Xn;
        end
        e0_sr(:, i) = e0_sr(:, i) + e0.^2;

        % --- Variant 1: mu(n) = alpha*log10(1 + (1/2)*(e(n)/delta)^2) ------
        % Depends on the current error only.
        for n = L:N-1
            Xn     = dlag(n:-1:n-L+1, i);
            y1(n)  = Xn' * W1(:, n);
            e1(n)  = d(n, i) - y1(n);
            mu1(n) = alpha * log10(1 + 0.5 * e1(n)^2 / delta^2);
            W1(:, n+1) = W1(:, n) + 2 * mu1(n) * e1(n) * Xn;
        end
        e1_sr(:, i) = e1_sr(:, i) + e1.^2;

        % --- Variant 2: mu(n) = alpha*log10(1 + (1/2)*|e(n)e(n-1)|/delta^2)
        % Uses the product of successive errors, so uncorrelated noise
        % tends to cancel and the step size reacts to the trend rather
        % than to a single noisy sample.
        for n = L:N-1
            Xn     = dlag(n:-1:n-L+1, i);
            y2(n)  = Xn' * W2(:, n);
            e2(n)  = d(n, i) - y2(n);
            mu2(n) = alpha * log10(1 + 0.5 * abs(e2(n)*e2(n-1)) / delta^2);
            W2(:, n+1) = W2(:, n) + 2 * mu2(n) * e2(n) * Xn;
        end
        e2_sr(:, i) = e2_sr(:, i) + e2.^2;
    end

    % --- Per-SNR diagnostics -----------------------------------------------
    figure;
    subplot(3,1,1); plot(nn, x);
    title('Input signal x(n)'); xlabel('Sample (n)'); ylabel('Amplitude');

    subplot(3,1,2); plot(nn, d(:,i));
    title(sprintf('Noisy signal d(n), SNR = %d dB', SNR(i)));
    xlabel('Sample (n)'); ylabel('Amplitude');

    subplot(3,1,3); plot(nn, y0, '-k', nn, y1, '-b', nn, y2, '-r');
    legend('LMS', 'Variant 1', 'Variant 2');
    title('Filter output y(n)'); xlabel('Sample (n)'); ylabel('Amplitude');

    X  = freqz(x, 1, N/2);
    D  = freqz(d(:,i), 1, N/2);
    Y0 = freqz(y0, 1, N/2);
    Y1 = freqz(y1, 1, N/2);
    [Y2, w] = freqz(y2, 1, N/2);

    figure;
    plot(w/pi, 20*log10(abs([X D Y0 Y1 Y2])));
    legend('x', 'd', 'LMS', 'Variant 1', 'Variant 2');
    title(sprintf('Spectra, SNR = %d dB', SNR(i)));
    xlabel('Normalised frequency (\times\pi rad/sample)'); ylabel('Magnitude (dB)');
end

% --- Averaged learning curves ---------------------------------------------
e0_sr = e0_sr / num_runs;
e1_sr = e1_sr / num_runs;
e2_sr = e2_sr / num_runs;

for i = 1:num_SNR
    figure;
    semilogy(nn, [e0_sr(:,i) e1_sr(:,i) e2_sr(:,i)]);
    legend('LMS', 'Variant 1', 'Variant 2');
    title(sprintf('Mean-square error, SNR = %d dB (%d runs)', SNR(i), num_runs));
    xlabel('Sample (n)'); ylabel('MSE');
end
