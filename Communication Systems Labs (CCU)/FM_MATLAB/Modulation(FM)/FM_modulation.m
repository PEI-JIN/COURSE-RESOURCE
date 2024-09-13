%1.讀取sample3.wav檔
[y,FS] = audioread('sample3.wav',[1, 20*8e3]);

%2.設定Fs=8k, FSC=64k, fc=16k
%3.將聲音做超取樣八倍為64k,並設定時間index t
Fs = 8e3;
Fsc = 8*8e3;
fc = 16e3;

%4.把聲音去除直流成分
z = resample(y,8,1);
t = 1:length(z);
z = z-sum(z)/length(z);

%事先宣告一個零矩陣，增進程式效能
q = zeros(length(z),1);

%6.積分m(t)
q(1) = z(1);
for i = 2: length(z)
    q(i) = q(i-1) + z(i);
end

%7.積分出來的向量/最大元素的絕對值 
q = q/max(abs(q));

%8.寫出FM調變公式
rxFM = cos(2*pi*fc*t'/Fsc + q(t));

%9.儲存FM訊號
save rxFM; 