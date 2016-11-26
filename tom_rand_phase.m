function [ image ] = tom_rand_phase( image )

phases = rand(size(image)).*pi.*2;
phases = complex(cos(phases), sin(phases));

imageft = fftn(image);
imageft = imageft.*phases;

image = real(ifftn(imageft));

end

