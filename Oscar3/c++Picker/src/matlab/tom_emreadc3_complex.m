function vol = tom_emreadc3_complex(name_re, name_im)

a_re = tom_emreadc3(name_re);
a_im = tom_emreadc3(name_im);
vol = a_re.Value + j*a_im.Value;
