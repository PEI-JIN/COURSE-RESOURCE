%%%�e�X�i��
load rxFM; %���JFM�T��
figure(1);
clf;
offset = 100000;
subplot(2,1,1);
plot(z(offset+1:offset+100)); %�e�X�T���T��
xlabel('Samples');
ylabel('Amplitude');

subplot(2,1,2);
plot(rxFM(offset+1:offset+100)); %�e�XFM�T��
xlabel('Samples');
ylabel('Amplitude');

%�]�w
%A.�C�q�o�i������(�o�̳]�w300)
%B.�C�q�o�i���I���W�v(�]�w5k)
%C.�T�������W�v(64k)
for i = 1:300
    a(i) = sinc(5e3 * (i-300/2)/64e3);
end

%�]�w�L�����ߪi�T�� 
B = [1,-1];

%�i��L��
y = conv(rxFM,B);

%���i��y
x = abs(y);

%�g�L�C�q�o�i��
result = conv(x,a);

%�N�n������X��
soundsc(result,64000);
