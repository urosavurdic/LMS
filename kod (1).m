close all
clear
clc
% Parametri simulacije
L = 100;     % Red filtra
N = 10000;   % Broj odbiraka
mi0=1e-5;
alpha =mi0;
delta = 0.01;
br_sim = 100;
kas=1;

% Ulazni signal X(n) = sin(2*pi*n/10)
nn = 0:N-1;
x = sin(2*pi*nn/20);

% Beli šum n_1(n) za različite SNR vrednosti (10dB, 20dB, 30dB)
SNR = [10, 20, 30];
br_SNR = length(SNR);
%n1 = zeros(N, br_SNR);

e0_sr = zeros(N, br_SNR);
e1_sr = zeros(N, br_SNR);
e2_sr = zeros(N, br_SNR);

for i = 1:br_SNR
    SNR_dB = SNR(i);
    d(:, i) = awgn(x, SNR_dB, 'measured');
    dkas(:, i)=[zeros(kas,1); d(1:end-kas,i)];

% Generisanje istog željenog signala d(n) = X(n) + n1(n) za sve SNR vrednosti
%d = X' + n1(:, 1);

% Inicijalizacija težinskih koeficijenata za oba adaptivna filtra



% Simulacija za oba adaptivna filtra
for iter = 1:br_sim
    W0 = zeros(L, N);
    W1 = zeros(L, N);
    W2 = zeros(L, N);
    % Prvi adaptivni filter sa μ(n) = alpha * log(1 + 1/2 * (e(n) * e(n-1)) / delta^2)
    y0 = zeros(N, 1);
    e0 = zeros(N, 1);
%    W1(:, iter+1) = W1(:, iter);
    for n = L:N-1
        Xn = dkas(n:-1:n-L+1,i);
        y0(n) = Xn' * W0(:, n);
        e0(n) = d(n) - y0(n); % Koristimo isti d_signal za sve SNR vrednosti
        W0(:, n+1) = W0(:, n) + 2 * mi0 * e0(n) * Xn;
    end
    e0_sr(:,i)=e0_sr(:,i)+e0.^2;

    y1 = zeros(N, 1);
    e1 = zeros(N, 1);
%    W1(:, iter+1) = W1(:, iter);
    for n = L:N-1
        Xn = dkas(n:-1:n-L+1,i);
        y1(n) = Xn' * W1(:, n);
        e1(n) = d(n) - y1(n); % Koristimo isti d_signal za sve SNR vrednosti
        mi1(n)=alpha*log10(1+0.5*e1(n)^2/delta^2);
        W1(:, n+1) = W1(:, n) + 2 * mi1(n) * e1(n) * Xn;
    end
    e1_sr(:,i)=e1_sr(:,i)+e1.^2;

    % Drugi adaptivni filter sa μ(n) = alpha * log(1 + 1/2 * (e(n) / delta)^2)
    y2 = zeros(N, 1);
    e2 = zeros(N, 1);
%    W2(:, iter+1) = W2(:, iter);
    for n = L:N-1
        Xn = dkas(n:-1:n-L+1,i);
        y2(n) = Xn' * W2(:, n);
        e2(n) = d(n) - y2(n); % Koristimo isti d_signal za sve SNR vrednosti
        mi2(n)=alpha*log10(1+0.5*abs(e2(n)*e2(n-1))/delta^2);
        W2(:, n+1) = W2(:, n) + 2 * mi2(n) * e2(n) * Xn;
    end
    e2_sr(:,i)=e2_sr(:,i)+e2.^2;
end
figure,plot(nn,W0);
figure,plot(nn,W1);
figure,plot(nn,W2);
% Prikazivanje grafika signala X(n), d(n) i y(n) za svaku simulaciju
figure;
subplot(3, 1, 1);
plot(nn, x);
title('Ulazni signal X(n)');
xlabel('Odbirci (n)');
ylabel('Amplituda');

subplot(3, 1, 2);
plot(nn, d(:,i));
title('Za[umljen signal d(n)');
xlabel('Odbirci (n)');
ylabel('Amplituda');

subplot(3, 1, 3);
plot(nn, y0, '-k', nn, y1, '-b', nn, y2, '-r');
legend('LMS','Adaptivni filtar 1', 'Adaptivni filtar 2');
title('Dobijeni signal y(n)');
xlabel('Odbirci (n)');
ylabel('Amplituda');

X=freqz(x,1,N/2);
D=freqz(d(:,i),1,N/2);
Y0=freqz(y0,1,N/2);
Y1=freqz(y1,1,N/2);
[Y2,wosa]=freqz(y2,1,N/2);
figure,plot(wosa/pi,20*log10(abs([X D Y0 Y1 Y2])));
legend('X','D','Y0','Y1','Y2');

end
e0_sr=e0_sr/br_sim;
e1_sr=e1_sr/br_sim;
e2_sr=e2_sr/br_sim;
figure,semilogy(nn,([e0_sr(:,1) e1_sr(:,1) e2_sr(:,1)]));
figure,semilogy(nn,([e0_sr(:,2) e1_sr(:,2) e2_sr(:,2)]));
figure,semilogy(nn,([e0_sr(:,3) e1_sr(:,3) e2_sr(:,3)]));

