function fo = f_WL(h, g, fs)
    fo = (fs/(2*pi)) * atan2(real(sqrt(imag(h).^2 - abs(g).^2)), real(h));
end
