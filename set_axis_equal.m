function  [lim_min, lim_max]= set_axis_equal(signal)
% signal is complex
    mi_re = min(real(signal(:))); ma_re = max(real(signal(:)));
    mi_im = min(imag(signal(:))); ma_im = max(imag(signal(:)));
    mi = min([mi_re, mi_im]); ma = max([ma_re, ma_im]);
    lim_min = mi - abs(mi/3);
    lim_max = ma + abs(ma/3);
end