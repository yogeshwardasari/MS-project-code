# MS-project-code

clc
clear all
test=load('-ascii','C:\Users\YOGESHWAR\Documents\MATLAB\PROJECT\samplerun.mat');
for i=1:15
    if i==1
        time=test(1,:);
    else
        channel(i-1,:)=test(i,:);
    end
end
figure(1)
hold on
%CHN=[AF3; F7; F3; FC5; T7; P7; O1; O2; P8; T8; FC6; F4; F8; AF4];
%[m,n]=size(CHN);

for i=1:14
    figure(1)
    plot(time,channel(i,:));
    legend('AF3','F7','F3','FC5','T7','P7','O1','O2','P8','T8','FC6','F4','F8','AF4')
    xlabel('Time')
    ylabel('Amplitude')
end

Fs=128;

figure(2)
hold on

for i=1:14
    FFT(i,:)=fft(channel(i,:));
    b=FFT(i,:);
    n=length(b);
    f=0:128/(n-1):128;
    figure(2)  
    plot(f(1:end/2),abs(FFT(i,1:end/2))/n) 
    c=abs(FFT(i,:));
    xlabel('Frequency(Hz)')
    ylabel('Amplitude')
end
hold off


figure(3)
hold on

for i=1:14
    [pxx(:,i),f]=pwelch(channel(i,:),[],[],[],Fs,'power'); 
    plot(f,(pxx(:,i)));
%     plot(f,pow2db(pxx(:,i)));
    peakchannel(i)=max(pow2db(pxx(:,i)));
    power(i)=sum((pxx(:,i)));
    powerdb(i)=(10*log10(power(i)));
    spec_en(i) = sum(pentropy(channel(i,:),Fs));
    frequency1=medfreq(pxx,f);
    frequency2=meanfreq(pxx,f);
    xlabel('Frequency(Hz)')
    ylabel('PSD (db/Hz)')
end

 k=1;
 m=1;
%% alpha power band
for i=1:14
    k=1;
     for j=1:length(pxx(:,i))
        if f(j)>=8 && f(j)<=13
            A_powerband(k,m)=pxx(j,i);
            g(k,1)=f(j);
            k=k+1; 
        end
     end
    A_BAND_power(i)=sum((A_powerband(:,i)));
    A_BAND_powerdb(i)=(10*log10(A_BAND_power(i)));
    frequency3=medfreq(A_powerband, g);
    frequency4=meanfreq(A_powerband, g);
m=m+1;
figure(4)
 plot(g,pow2db(A_powerband))
 peakchannel2(i)=max(pow2db(A_BAND_power));
 xlabel('Frequency(Hz)')
 ylabel('Alpha PSD(db/Hz)')
end

k1=1;
 m1=1;
%% Theta power band
for i=1:14
    k1=1;
     for j=1:length(pxx(:,i))
        if f(j)>=4 && f(j)<=8
            A1_powerband(k1,m1)=pxx(j,i);
            g1(k1,1)=f(j);
            k1=k1+1; 
        end
     end
    A1_BAND_power(i)=sum((A1_powerband(:,i)));
    A1_BAND_powerdb(i)=(10*log10(A1_BAND_power(i)));
    frequency5=medfreq(A1_powerband, g1);
    frequency6=meanfreq(A1_powerband, g1);
m1=m1+1;
figure(5)
 plot(g1,pow2db(A1_powerband))
 peakchannel3(i)=max(pow2db(A1_BAND_power));
 xlabel('Frequency(Hz)')
 ylabel('Theta PSD(db/Hz)')
end
