%1.Ū��sample3.wav��
[y,FS] = audioread('sample3.wav',[1, 20*8e3]);

%2.�]�wFs=8k, FSC=64k, fc=16k
%3.�N�n�����W���ˤK����64k,�ó]�w�ɶ�index t
Fs = 8e3;
Fsc = 8*8e3;
fc = 16e3;

%4.���n���h�����y����
z = resample(y,8,1);
t = 1:length(z);
z = z-sum(z)/length(z);

%�ƥ��ŧi�@�ӹs�x�}�A�W�i�{���į�
q = zeros(length(z),1);

%6.�n��m(t)
q(1) = z(1);
for i = 2: length(z)
    q(i) = q(i-1) + z(i);
end

%7.�n���X�Ӫ��V�q/�̤j����������� 
q = q/max(abs(q));

%8.�g�XFM���ܤ���
rxFM = cos(2*pi*fc*t'/Fsc + q(t));

%9.�x�sFM�T��
save rxFM; 