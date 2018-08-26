function J = rk4exp2(w,dt,g,E,KT)
E2=E.^2;
a=g.*reshape((fft2(real(ifft2(reshape(w.',KT,KT).')).^2))',KT^2,1); %you have to reshape w to ifft2
b=g.*reshape((fft2(real(ifft2(reshape((E.*(w+a/2)).',KT,KT).')).^2))',KT^2,1);
c=g.*reshape((fft2(real(ifft2(reshape((E.*w+b/2).',KT,KT).')).^2))',KT^2,1);
d=g.*reshape((fft2(real(ifft2(reshape((E2.*w+E.*c).',KT,KT).')).^2))',KT^2,1);

J=E2.*w+(E2.*a+2*E.*(b+c)+d)/6;
end