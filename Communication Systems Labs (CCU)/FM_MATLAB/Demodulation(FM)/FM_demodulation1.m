%%%畫出波形
load rxFM; %載入FM訊號
figure(1);
clf;
offset = 100000;
subplot(2,1,1);
plot(z(offset+1:offset+100)); %畫出訊息訊號
xlabel('Samples');
ylabel('Amplitude');

subplot(2,1,2);
plot(rxFM(offset+1:offset+100)); %畫出FM訊號
xlabel('Samples');
ylabel('Amplitude');

%設定
%A.低通濾波器長度(這裡設定300)
%B.低通濾波器截止頻率(設定5k)
%C.訊號取樣頻率(64k)
for i = 1:300
    a(i) = sinc(5e3 * (i-300/2)/64e3);
end

%設定微分器脈波響應 
B = [1,-1];

%進行微分
y = conv(rxFM,B);

%全波整流
x = abs(y);

%經過低通濾波器
result = conv(x,a);

%將聲音播放出來
soundsc(result,64000);
