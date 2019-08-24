%% Mid-Term Project(Extra Credit)
%% ALC 5564/5764
% Mohak Bheda
%% Section 1

%To construct excitation for System ID
%To compute and plot properly scaled power spectrum for excitation
%Plot excitation in time domain to show that it does not saturate the DAC

%To Display: 1 Figure only

clc;
clear all;

fs = 20;    %Sample rate
ts = 1/fs;
N = 72000;  %Samples %Reduced samples-Extra credit
wn= [0.01 9.9];
t = (0:(N-1))*ts;

scale = 2.5;
[num,den] = butter(4,wn/(fs/2));    %Filter of Order 4
u_temp = scale * filter(num,den,randn(N,2));
u = u_temp';
u1 = u(1,:);
u2 = u(2,:);

nfft = 2^11;
wndo = nfft;
ovlp = nfft/2;

[Suu1,fr] = cpsd(u1,u1,wndo,ovlp,nfft,fs);
[Suu2,fr] = cpsd(u2,u2,wndo,ovlp,nfft,fs);

figure()
semilogx(fr, 10 *log10(abs(Suu1)*fs/2),'lineWidth',1.5);
hold on
semilogx(fr, 10 *log10(abs(Suu2)*fs/2),'lineWidth',1.5);
hold off
legend('Input - U1', 'Input - U2','location','southeast');
title('Power Spectrum for Excitation')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
xlim([0.01 10])

figure()
plot(t,u)
legend('Input - U1', 'Input - U2');
title('Excitation in time domain')
xlabel('time')
ylabel('Excitation')
xlim([0 N/fs])
%% *Section 2*
%%
%Apply excitations to the dynamic system to obtain noisy response signals
%Estimate SNR for each of the two paths
%Note1:Better System ID results will be obtained if SNR is maximized
%Note2:Do not saturate the excitations

%To Display: Properly scaled SNR results 

y= s19_plant(u);
y1 = s19_plant([u(1,:);zeros(1,N)]); %Output 1, with only input 1
y2 = s19_plant([zeros(1,N);u(2,:)]); %Output 2, with only input 2

noise_input = zeros(2,N);   %Noise input 
n = s19_plant(noise_input);
n1 = n(1,:);
n2 = n(2,:);

sigma_y11 = std(y1(1,:));
sigma_y21 = std(y2(1,:));
sigma_y12 = std(y1(2,:));
sigma_y22 = std(y2(2,:));

sigma_n1 = std(n1);
sigma_n2 = std(n2);


SNR1 = 10 * log10((sigma_y11^2-sigma_n1^2)/sigma_n1^2); %SNR for O/P 1, I/P 1
SNR2 = 10 * log10((sigma_y21^2-sigma_n1^2)/sigma_n1^2); %SNR for O/P 2, I/P 1
SNR3 = 10 * log10((sigma_y12^2-sigma_n2^2)/sigma_n2^2); %SNR for O/P 1, I/P 2
SNR4 = 10 * log10((sigma_y22^2-sigma_n2^2)/sigma_n2^2); %SNR for O/P 2, I/P 2

display(SNR1)
display(SNR2)
display(SNR3)
display(SNR4)
%% *Section 3*
%%
%To plot properly scaled power spectrum for the responses
%On the same plot, plot power spectrum for the noise signals
%Plot responses in time domain to show that they do not saturate the ADC's

%To Display: 4 separate axes(2 for freq, 2 for time) in a single 2X2 Figure


[Syy1,fr] = cpsd(y(1,:),y(1,:),wndo,ovlp,nfft,fs);  %Power Spectrum for response 1
[Syy2,fr] = cpsd(y(2,:),y(2,:),wndo,ovlp,nfft,fs);  %Power Spectrum for response 2
[Snn1,fr] = cpsd(n1,n1,wndo,ovlp,nfft,fs);          %Power Spectrum for noise 1
[Snn2,fr] = cpsd(n2,n2,wndo,ovlp,nfft,fs);          %Power Spectrum for noise 2


figure()

subplot(2,2,1);
semilogx(fr, 10 *log10(abs(Syy1)*fs/2),'lineWidth',1.5);
hold on
semilogx(fr, 10 *log10(abs(Snn1)*fs/2),'lineWidth',1.5);
hold off
grid on;
legend('Output - y1', 'Noise - n1','location','southeast');
set(gcf, 'units', 'points', 'position',[0,0,500,400])
ylim([-100 30])
title('Power Spectrum- y1,n1')
xlabel('Frequency')
ylabel('Magnitude')


subplot(2,2,2);
semilogx(fr, 10 *log10(abs(Syy2)*fs/2),'lineWidth',1.5);
hold on
semilogx(fr, 10 *log10(abs(Snn2)*fs/2),'lineWidth',1.5);
hold off
grid on;
legend('Output - y2', 'Noise - n2','location','southeast');
set(gcf, 'units', 'points', 'position',[0,0,500,400])
ylim([-100 30])
title('Power Spectrum- y2,n2')
xlabel('Frequency')
ylabel('Magnitude')

subplot(2,2,3);
plot(t, y(1,:))
hold on
plot(t, n1,'lineWidth',1.5)
hold off
grid on;
legend('Output - y1', 'Noise - n1','location','southeast');
set(gcf, 'units', 'points', 'position',[0,0,500,400])
title('Power Spectrum time domain- y1,n1')
xlabel('Frequency')
ylabel('Magnitude')

subplot(2,2,4);
plot(t, y(2,:))
hold on
plot(t, n2,'lineWidth',1.5)
hold off
grid on;
ylim([-2,2])
legend('Output - y2', 'Noise - n2','location','southeast');
set(gcf, 'units', 'points', 'position',[0,0,500,400])
title('Power Spectrum time domain- y2,n2')
xlabel('Frequency')
ylabel('Magnitude')
%% *Section 4*
%%
%To apply an apt Freq Domain Estimation Technique to estimate Freqresponse
%Do this for each open loop path
%Also, estimate coherence for each path

%To Display: Nothing for this section

[sys_y11_est,fr] = tfestimate(u1,y1(1,:),wndo,ovlp,nfft,fs); %Estimation for Input 1 Output 1
[sys_y21_est,fr] = tfestimate(u1,y1(2,:),wndo,ovlp,nfft,fs); %Estimation for Input 1 Output 2
[sys_y12_est,fr] = tfestimate(u2,y2(1,:),wndo,ovlp,nfft,fs); %Estimation for Input 2 Output 1
[sys_y22_est,fr] = tfestimate(u2,y2(2,:),wndo,ovlp,nfft,fs); %Estimation for Input 2 Output 2

gamma11 = mscohere(u1,y1(1,:),wndo, ovlp, nfft, fs);  %Coherence for Input 1 Output 1
gamma21 = mscohere(u1,y1(2,:),wndo, ovlp, nfft, fs);  %Coherence for Input 1 Output 2
gamma12 = mscohere(u2,y2(1,:),wndo, ovlp, nfft, fs);  %Coherence for Input 2 Output 1
gamma22 = mscohere(u2,y2(2,:),wndo, ovlp, nfft, fs);  %Coherence for Input 2 Output 2
%% Section 5
%%
%First, Determine Disc Time TF models for the two paths; Use 'invfreqz' Func
%Choose same No. of Poles for each path
%Store Num & Den polynomials in cell arrays with Apt dimensions
%Convert Num & Den cell arrays into a Disc Time TF LTI object
%Convert the LTI object into a minimum realization; Use 'minreal' Func

%To Display: Nothing for this section

[Suu1,fr] = cpsd(u1,u1,wndo,ovlp,nfft,fs);
[Suu2,fr] = cpsd(u2,u2,wndo,ovlp,nfft,fs);

Wt11 = Suu1;
Wt21 = Suu2;
Wt12 = Suu1;
Wt22 = Suu2;

Nb = 17; %Poles
Na = 17; %Zeros

[B11,A11] = invfreqz(sys_y11_est, 2*pi*fr/fs, Nb, Na, Wt11);
[B21,A21] = invfreqz(sys_y21_est, 2*pi*fr/fs, Nb, Na, Wt21);
[B12,A12] = invfreqz(sys_y12_est, 2*pi*fr/fs, Nb, Na, Wt12);
[B22,A22] = invfreqz(sys_y22_est, 2*pi*fr/fs, Nb, Na, Wt22);

Num_Cell = {B11, B21, B12, B22}
Den_Cell = {A11, A21, A12, A22}

sys = tf(Num_Cell, Den_Cell, ts); %System trasfer function in 1X4 block
sys_y11 = sys(1,1,:); %System TF for input 1 output 1
sys_y21 = sys(1,2,:); %System TF for input 1 output 2
sys_y12 = sys(1,3,:); %System TF for input 2 output 1
sys_y22 = sys(1,4,:); %System TF for input 2 output 2

sys_mr = minreal(sys); %minimum realization for the system
sys_mr_y11 = minreal(sys_y11); %minimum realization for input 1 output 1
sys_mr_y21 = minreal(sys_y21); %minimum realization for input 1 output 2
sys_mr_y12 = minreal(sys_y12); %minimum realization for input 2 output 1
sys_mr_y22 = minreal(sys_y22); %minimum realization for input 2 output 2
%% Section 6
%%
%First, Convert the min realization TF LTI object to a disc Time SS LTI object
%I/P the SS model to the 'balreal' Func to generate a balanced realization 
%Extract the Hankel singular values form 'balreal' to find Apt cutoff thres
%Generate a reduced order SS Disc Time LTI model; Use 'modred' Func

%To Display: 1 figure only-Bar chart(For HSV)

sys_ss = ss(sys_mr); %State space model
[sysb, g] = balreal(sys_ss); %balanced realization for the SS model
figure()
%hsvd(sysb,opts);

bar(g)
set(gca, 'YScale', 'log');
elim = (g < 1 | isinf(g));
xlabel('X axis')
ylabel('Y axis')

sys_red = modred(sysb, elim, 'del'); %reduced order of the system
sys_red_y11 = sys_red(1,1,:);
sys_red_y21 = sys_red(1,2,:);
sys_red_y12 = sys_red(1,3,:);
sys_red_y22 = sys_red(1,4,:);
%% Section 7
%%
%First, Compute z-domain Eigenvalues(poles) of the LTI object from Sec 6
%For Every Eigenvalue, Compute natural Freqs(in Hz)
%For Every Eigenvalue, Compute Damping ratios-damp func

%To Display: Numerical values

Pole_red = pole(sys_red); %Eigenvalue of the LTI object
[wn_red,zeta_red, pole_red] = damp(sys_red);
damp(sys_red) %Computing natural frequencies and Damping
%% Section 8
%%
%First, Compute Eigenvalues for open loop TF model in Section 5
%Then using 'zgrid', generate a z-domain grid and plot each set of poles 
%Three set of poles:- For Y1 path, Y2 path, and for Disc time SS model(Sec 6)
%Note: Use clearly labelled markers for each path + Use 'axis equal'

%To Display: 1 Figure Only- Plot for each set of Poles

Pole_mr_y11  =  pole(sys_mr_y11); %Pole for min realization for I/P 1, O/P 1
Pole_mr_y21  =  pole(sys_mr_y21); %Pole for min realization for I/P 1, O/P 2

Pole_mr_y12  =  pole(sys_mr_y12); %Pole for min realization for I/P 2, O/P 1
Pole_mr_y22  =  pole(sys_mr_y22); %Pole for min realization for I/P 2, O/P 2


Pole_sys_red_y11 = pole(sys_red_y11); %Pole for reduced order for I/P 1, O/P 1
Pole_sys_red_y21 = pole(sys_red_y21); %Pole for reduced order for I/P 1, O/P 2
Pole_sys_red_y12 = pole(sys_red_y12); %Pole for reduced order for I/P 2, O/P 1
Pole_sys_red_y22 = pole(sys_red_y22); %Pole for reduced order for I/P 2, O/P 2

Pole_sys_red = pole(sys_red); 

figure()

plot(Pole_mr_y11, 'o');
hold on
plot(Pole_mr_y21, 'x');
plot(Pole_mr_y12, '+');
plot(Pole_mr_y22, '*');
plot(Pole_sys_red_y11,'o');
plot(Pole_sys_red_y21,'x');
plot(Pole_sys_red_y12,'+');
plot(Pole_sys_red_y22,'*');
zgrid;
legend('unit circle','zeta','wn','y11','y21','y12','y22','location','bestoutside')
xlabel('Real axis')
ylabel('Imaginary Axis')
axis equal;
%% Section 9
%%
%Compute Freq Response of LTI model from Sec 5; Use 'freqresp' Func
%Compute Freq Response of LTI model from Sec 6; Use 'freqresp' Func

%To Display: Nothing for this section


FRF_mr = squeeze(freqresp(sys_mr,2*pi*fr)); %Freq response for min real model

FRF_red_y11 = squeeze(freqresp(sys_red_y11,2*pi*fr)); %Freq response for reduced order model
FRF_red_y21 = squeeze(freqresp(sys_red_y21,2*pi*fr)); 
FRF_red_y12 = squeeze(freqresp(sys_red_y12,2*pi*fr));
FRF_red_y22 = squeeze(freqresp(sys_red_y22,2*pi*fr));
%% Section 10
%%
%Generate a dB mag vs. log Freq Plot for each path
%Generate a Phase(degrees) vs. log Freq Plot for each path
%Also plot the coherence(Sec 4) for each path

%To display: 2 figures with 3 axes in a 3X1 grid;must be properly labelled

%fr = logspace(-2,2,1000); % frequency vector (Hz)


%figure()


figure;

subplot(3,1,1)
semilogx(fr,20*log10(abs(sys_y11_est)),'lineWidth',1.5);
hold on;
semilogx(fr,20*log10(abs(squeeze(freqresp(sys_y11,2*pi*fr)))),'lineWidth',1.5);
semilogx(fr,20*log10(abs(squeeze(freqresp(sys_red_y11,2*pi*fr)))),'lineWidth',1.5);
grid on;
legend('TF', 'TF-Min real','TF-Reduced order','location','southwest');
title('Mag vs Freq plot for Y1,U1');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3,1,2)
semilogx(fr,(180/pi)*angle(sys_y11_est),'lineWidth',1.5);
hold on;
semilogx(fr,(180/pi)*angle(squeeze(freqresp(sys_y11,2*pi*fr))),'lineWidth',1.5);
semilogx(fr,(180/pi)*angle(squeeze(freqresp(sys_red_y11,2*pi*fr))),'lineWidth',1.5);
grid on;
legend('TF', 'TF-Min real','TF-Reduced order','location','northwest');
title('Phase(deg) vs Freq plot for Y1,U1');
xlabel('Frequency (Hz)');
ylabel('Phase');
ylim([-180 180])

subplot(3,1,3)
semilogx(fr, gamma11,'lineWidth',1.5);
title('Coherence vs Freq plot for Y1,U1');
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;

figure;

subplot(3,1,1)
semilogx(fr,20*log10(abs(sys_y21_est)),'lineWidth',1.5);
hold on;
semilogx(fr,20*log10(abs(squeeze(freqresp(sys_y21,2*pi*fr)))),'lineWidth',1.5);
semilogx(fr,20*log10(abs(squeeze(freqresp(sys_red_y21,2*pi*fr)))),'lineWidth',1.5);
grid on;
legend('TF', 'TF-Min real','TF-Reduced order','location','southwest');
title('Mag vs Freq plot for Y2,U1');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3,1,2)
semilogx(fr,(180/pi)*angle(sys_y21_est),'lineWidth',1.5);
hold on;
semilogx(fr,(180/pi)*angle(squeeze(freqresp(sys_y21,2*pi*fr))),'lineWidth',1.5);
semilogx(fr,(180/pi)*angle(squeeze(freqresp(sys_red_y21,2*pi*fr))),'lineWidth',1.5);
grid on;
legend('TF', 'TF-Min real','TF-Reduced order','location','northwest');
title('Phase(deg) vs Freq plot for Y2,U1');
xlabel('Frequency (Hz)');
ylabel('Phase');
ylim([-180 180])

subplot(3,1,3)
semilogx(fr, gamma21 ,'lineWidth',1.5);
title('Coherence vs Freq plot for Y2,U1');
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;

figure;

subplot(3,1,1)

semilogx(fr,20*log10(abs(sys_y12_est)),'lineWidth',1.5);
hold on;
semilogx(fr,20*log10(abs(squeeze(freqresp(sys_y12,2*pi*fr)))),'lineWidth',1.5);
semilogx(fr,20*log10(abs(squeeze(freqresp(sys_red_y12,2*pi*fr)))),'lineWidth',1.5);
grid on;
legend('TF', 'TF-Min real','TF-Reduced order','location','southwest');
title('Mag vs Freq plot for Y1,U2');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3,1,2)
semilogx(fr,(180/pi)*angle(sys_y12_est),'lineWidth',1.5);
hold on;
semilogx(fr,(180/pi)*angle(squeeze(freqresp(sys_y12,2*pi*fr))),'lineWidth',1.5);
semilogx(fr,(180/pi)*angle(squeeze(freqresp(sys_red_y12,2*pi*fr))),'lineWidth',1.5);
grid on;
legend('TF', 'TF-Min real','TF-Reduced order','location','southwest');
title('Phase(deg) vs Freq plot for Y1,U2');
xlabel('Frequency (Hz)');
ylabel('Phase');
ylim([-180 180])

subplot(3,1,3)
semilogx(fr, gamma12,'lineWidth',1.5);
title('Coherence vs Freq plot for Y1,U2');
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;

figure;
subplot(3,1,1)
semilogx(fr,20*log10(abs(sys_y22_est)),'lineWidth',1.5);
hold on;
semilogx(fr,20*log10(abs(squeeze(freqresp(sys_y22,2*pi*fr)))),'lineWidth',1.5);
semilogx(fr,20*log10(abs(squeeze(freqresp(sys_red_y22,2*pi*fr)))),'lineWidth',1.5);
grid on;
legend('TF', 'TF-Min real','TF-Reduced order','location','southwest');
title('Mag vs Freq plot for Y2,U2');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3,1,2)
semilogx(fr,(180/pi)*angle(sys_y22_est),'lineWidth',1.5);
hold on;
semilogx(fr,(180/pi)*angle(squeeze(freqresp(sys_y22,2*pi*fr))),'lineWidth',1.5);
semilogx(fr,(180/pi)*angle(squeeze(freqresp(sys_red_y22,2*pi*fr))),'lineWidth',1.5);
grid on;
legend('TF', 'TF-Min real','TF-Reduced order','location','northwest');
title('Phase(deg) vs Freq plot for Y2,U2');
xlabel('Frequency (Hz)');
ylabel('Phase');
ylim([-180 180])

subplot(3,1,3)
semilogx(fr, gamma22,'lineWidth',1.5);
title('Coherence vs Freq plot for Y2,U2');
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;