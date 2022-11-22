function fo = f_SL(h, fs)
    fo = (fs/(2*pi)) * atan2(imag(h), real(h));
end