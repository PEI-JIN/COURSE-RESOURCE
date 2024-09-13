clc;
clear all;
load rxQAM16n;  %載入QAM訊號

%%計算alpha部分%%
alpha = (sum(rxQAM16(1:32).*phi2)+sum(rxQAM16(1:32).*phi1))/2;

%%星象圖部分%%
for i=1:length(rxQAM16)/32-1
    c1(i)=sum(rxQAM16(i*32+1:i*32+32).*phi2)/alpha; %AR
    c2(i)=sum(rxQAM16(i*32+1:i*32+32).*phi1)/alpha; %AI
    
    if c1(i) < -2                        %AR
        m((i-1)*4+1:(i-1)*4+2)=[0 0];
    elseif c1(i) < 0
        m((i-1)*4+1:(i-1)*4+2)=[0 1];
    elseif c1(i) < 2
        m((i-1)*4+1:(i-1)*4+2)=[1 0];
    else
        m((i-1)*4+1:(i-1)*4+2)=[1 1];
    end
    
    if c2(i) < -2                        %AI
        m((i-1)*4+3:(i-1)*4+4)=[0 0];
    elseif c2(i) < 0
        m((i-1)*4+3:(i-1)*4+4)=[0 1];
    elseif c2(i) < 2
        m((i-1)*4+3:(i-1)*4+4)=[1 0];
    else
        m((i-1)*4+3:(i-1)*4+4)=[1 1];
    end
end

for i=1:length(rxQAM16N)/32-1
    c1(i)=sum(rxQAM16N(i*32+1:i*32+32).*phi2)/alpha; %AR
    c2(i)=sum(rxQAM16N(i*32+1:i*32+32).*phi1)/alpha; %AI
    
    if c1(i) < -2                        %AR
        m((i-1)*4+1:(i-1)*4+2)=[0 0];
    elseif c1(i) < 0
        m((i-1)*4+1:(i-1)*4+2)=[0 1];
    elseif c1(i) < 2
        m((i-1)*4+1:(i-1)*4+2)=[1 0];
    else
        m((i-1)*4+1:(i-1)*4+2)=[1 1];
    end
    
    if c2(i) < -2                        %AI
        m((i-1)*4+3:(i-1)*4+4)=[0 0];
    elseif c2(i) < 0
        m((i-1)*4+3:(i-1)*4+4)=[0 1];
    elseif c2(i) < 2
        m((i-1)*4+3:(i-1)*4+4)=[1 0];
    else
        m((i-1)*4+3:(i-1)*4+4)=[1 1];
    end
end

plot(c1,c2,'*');
grid on

mf = zeros(1,length(m)/8);

for n = 1:1:length(m)/8
   mf(n) =  m(8*n)*16 + m(8*n-1)*32 +m(8*n-2)*64 + m(8*n-3)*128 + m(8*n-4)*1 + m(8*n-5)*2 + m(8*n-6)*4 + m(8*n-7)*8;
end

%播放聲音
soundsc(mf,8e3);