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
%A.�C�q���i������(�o�̳]�w300)
%B.�T�������W�v(64k)

%�Ĥ@�C�q�o�i�I���W�v20K
for i = 1:300
    a1(i) = sinc(20e3 * (i-300/2)/64e3);
end

%�ĤG�C�q�o�i�I���W�v5K
for i = 1:300
    a2(i) = sinc(5e3 * (i-300/2)/64e3);
end

%�NFM�T�����O���W cos(2pi*fc*t/�����W�v) �� sin(2pi*fc*t /�����W�v)
b1 = rxFM.*cos(2*pi*fc*t'/64e3);
b2 = rxFM.*sin(2*pi*fc*t'/64e3);

%���O�g�L�C�q�o�i��1
x1 = conv(a1,b1);
x2 = conv(a1,b2);

%�]�@�Ƽ� �N�곡���P�곡���۬ۥ[
y = x1 + j*x2;

%�ƼưT�����e�@�I(n)���W��@�I(n-1)���@�m
z = y(1:length(y)-1) .* conj( y(2:length(y)) );

%���n���q�ƼƬۦ줤���X
c = angle(z);

%����X���T���g�L�C�q�o�i��2
result = conv(c,a2);

%�����n��
soundsc(result,64000);
