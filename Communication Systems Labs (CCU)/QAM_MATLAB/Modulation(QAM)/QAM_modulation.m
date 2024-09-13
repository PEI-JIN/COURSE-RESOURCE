clear all;
am=[
-3,
-1,
1,
3
];

fc = 8e3*8; %¸üªi
fsample = 64e3*8; %¨ú¼ËÀW²v
brate = 2e3*8;
T = 1/brate;
Tsample = 1/fsample;
phi1 = sin(2*pi*fc*[0:31]*Tsample);
phi2 = cos(2*pi*fc*[0:31]*Tsample);
[y,FS]=audioread('sample.wav',[1+8000, 10*8e3+8000]);
% soundsc(y,FS,NBITS);

mean_y = mean(y);
max_y = max(abs(y-mean_y));
z = round((y-mean_y)/max_y*120+127);
rxQAM16=zeros(1, length (z)*2*32+32);
rxQAM16(1:32) = 1*phi1+1*phi2;

for i = 1: length (z)
rxQAM16((i-1)*64+1+32: (i-1)*64+32+32) = am(bitget(z(i),1)+2*bitget(z(i),2)+1)*phi1 + am(bitget(z(i),3) + 2*bitget(z(i),4)+1)*phi2;
rxQAM16((i-1)*64+33+32: (i-1)*64+64+32) = am(bitget(z(i),5)+ 2*bitget(z(i),6)+1)*phi1 + am(bitget(z(i),7) + 2*bitget(z(i),8)+1)*phi2;
end


rxQAM16N = rxQAM16+randn(size(rxQAM16))*0.4;
save rxQAM16N;