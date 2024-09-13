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
%A.低通路波器長度(這裡設定300)
%B.訊號取樣頻率(64k)

%第一低通濾波截止頻率20K
for i = 1:300
    a1(i) = sinc(20e3 * (i-300/2)/64e3);
end

%第二低通濾波截止頻率5K
for i = 1:300
    a2(i) = sinc(5e3 * (i-300/2)/64e3);
end

%將FM訊號分別乘上 cos(2pi*fc*t/取樣頻率) 及 sin(2pi*fc*t /取樣頻率)
b1 = rxFM.*cos(2*pi*fc*t'/64e3);
b2 = rxFM.*sin(2*pi*fc*t'/64e3);

%分別經過低通濾波器1
x1 = conv(a1,b1);
x2 = conv(a1,b2);

%設一複數 將實部項與虛部項相相加
y = x1 + j*x2;

%複數訊號的前一點(n)乘上後一點(n-1)的共軛
z = y(1:length(y)-1) .* conj( y(2:length(y)) );

%把聲音從複數相位中取出
c = angle(z);

%把取出的訊號經過低通濾波器2
result = conv(c,a2);

%播放聲音
soundsc(result,64000);
