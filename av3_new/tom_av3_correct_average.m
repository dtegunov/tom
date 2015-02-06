function avg_corr = tom_av3_correct_average(avg, compound_wedge)
    lowpass = floor(size(avg,1)/2)-2;
    compound_wedge = 1./compound_wedge;
    h = find(compound_wedge > 100000);
    compound_wedge(h) = 0;% take care for inf
    avg_corr = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(avg)).*compound_wedge,lowpass))));
end
