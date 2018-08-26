function J = rk4exp(w,dt,g,E2)
t=n*dt;%g=-.5i*dt*k;%kdv
%E=exp(dt*Lop/2); 
E2=E.^2;
a=g.*fft(real(ifft(w)).^2);
b=g.*fft(real(ifft(E.*(w+a/2))).^2);
c=g.*fft(real(ifft(E.*w+b/2)).^2);
d=g.*fft(real(ifft(E2.*w+E.*c)).^2);
J=E2.*w+(E2.*a+2*E.*(b+c)+d)/6;
end