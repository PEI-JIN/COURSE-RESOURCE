clc;
clear all;
load rxQAM16n;  %載入QAM訊號

rx11 = rxQAM16(1:32).*phi1;  %phi1 = sin(2*pi*fc*[0:31]*Tsample);
rx12 = rxQAM16(1:32).*phi2;  %phi2 = cos(2*pi*fc*[0:31]*Tsample);

q1 = zeros(length(rx11),1);
q1(1) = rx11(1); 
for i = 2:32
    q1(i) = q1(i-1) + rx11(i);
end

q2 = zeros(length(rx12),1);
q2(1) = rx12(1); 
for i = 2:32
    q2(i) = q2(i-1) + rx12(i);
end

alpha = (q1(32) + q2(32))/2;

A_R1 = zeros(1,length(z));
A_I1 = zeros(1,length(z));

A_R2 = zeros(1,length(z));
A_I2 = zeros(1,length(z));

A_R = zeros(1,2*length(z));
A_I = zeros(1,2*length(z));

for i = 1: length(z)
    rx1 = rxQAM16N((i-1)*64+1+32:(i-1)*64+32+32).*phi1;
    rx2 = rxQAM16N((i-1)*64+1+32:(i-1)*64+32+32).*phi2;
  
    q1 = zeros(length(rx1),1);
    q1(1) = rx1(1); 
    for j = 2:32
        q1(j) = q1(j-1) + rx1(j);
    end

    q2 = zeros(length(rx2),1);
    q2(1) = rx2(1); 
    for k = 2:32
        q2(k) = q2(k-1) + rx2(k);
    end
    
    A_R1(1,i) = q1(32)/alpha; %AR_1
    A_I1(1,i) = q2(32)/alpha; %AI_1
end

for i = 1: length(z)
    rx21 = rxQAM16N((i-1)*64+33+32:(i-1)*64+64+32).*phi1;
    rx22 = rxQAM16N((i-1)*64+33+32:(i-1)*64+64+32).*phi2;
  
    q1 = zeros(length(rx21),1);
    q1(1) = rx21(1); 
    for j = 2:32
        q1(j) = q1(j-1) + rx21(j);
    end

    q2 = zeros(length(rx22),1);
    q2(1) = rx22(1); 
    for k = 2:32
        q2(k) = q2(k-1) + rx22(k);
    end
    
    A_R2(1,i) = q1(32)/alpha; %AR_2
    A_I2(1,i) = q2(32)/alpha; %AI_2
end

for i = 1: 2*length(z)
    if i <= length(z)
        A_R(1,i) = A_R1(1,i); %AR_1
        A_I(1,i) = A_I1(1,i); %AI_1
    else
        A_R(1,i) = A_R2(1,i-80000); %AR_2
        A_I(1,i) = A_I2(1,i-80000); %AI_2
    end
end

%%%畫出 16-QAM Constellation
figure(1);
clf;
plot(A_I, A_R,'x');
axis([-4, 4, -4, 4]);
grid on

xlabel('A_R');
ylabel('A_I');
title('16-QAM constellation');

m = zeros(length(z),1);

for i = 1:length(z)
    temp = round(A_R1(1,i));
	switch temp 
        case -3
            m(i,1) = m(i,1);
        case -1
            m(i,1) = bitset(m(i,1),1);
        case 1
            m(i,1) = bitset(m(i,1),2);
        case 3
            m(i,1) = bitset(m(i,1),1);
            m(i,1) = bitset(m(i,1),2);
    end
    
    temp2 = round(A_I1(1,i));
    switch temp2 
        case -3
            m(i,1) = m(i,1);
        case -1
            m(i,1) = bitset(m(i,1),3);
        case 1
            m(i,1) = bitset(m(i,1),4);
        case 3
            m(i,1) = bitset(m(i,1),3);
            m(i,1) = bitset(m(i,1),4);
    end
    
    temp3 = round(A_R2(1,i));
    switch temp3
        case -3
            m(i,1) = m(i,1);
        case -1
            m(i,1) = bitset(m(i,1),5);
        case 1
            m(i,1) = bitset(m(i,1),6);
        case 3
            m(i,1) = bitset(m(i,1),5);
            m(i,1) = bitset(m(i,1),6);
    end
    
    temp4 =round(A_I2(1,i));
    switch temp4
        case -3
            m(i,1) = m(i,1);
        case -1
            m(i,1) = bitset(m(i,1),7);
        case 1
            m(i,1) = bitset(m(i,1),8);
        case 3
            m(i,1) = bitset(m(i,1),7);
            m(i,1) = bitset(m(i,1),8);
    end
end

%result = (m-127)/120*max_y + mean_y;

%播放聲音
soundsc(m,FS);
