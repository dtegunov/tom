function HdArray=MakeImagicHeader(n1,n2,nim)
% HdArray=MakeImagicHeader(n1,n2,nim)
% Create a default HdArray for a stack of nim
% n1 x n2 pixel images.  We assume the data will be float32,
% and will be written little-ended.
% We fill in the minimal information (image no., size, date)
% The HdArray is a struct with the following fields:
%   Vals [double(199,nim)]   Numeric values for header
%   Name [char(nim,80)]      Name strings (set to null here)
%   Strings [char(nim,228)]  Strings (set to null here)
%   Type='REAL'              Always set here
%   ByteOrder='ieee-le'      Always set here

HdArray.Vals=zeros(199,nim);
HdArray.Name=char(zeros(nim,80));
HdArray.Strings=char(zeros(nim,228));
HdArray.Type='REAL';
HdArray.ByteOrder='ieee-le';

t=clock;
HdArray.Vals(1,:)=1:nim; % image index
HdArray.Vals(2,:)=nim-(1:nim);  % number of images following

HdArray.Vals(4,:)=1;     % number of headers per image
HdArray.Vals(5,:)=t(3);  % day
HdArray.Vals(6,:)=t(2);  % month
HdArray.Vals(7,:)=t(1);  % year
HdArray.Vals(8,:)=t(4);  % hour
HdArray.Vals(9,:)=t(5);  % minute
HdArray.Vals(10,:)=round(t(6)); % second
HdArray.Vals(12,:)=n1*n2;
HdArray.Vals(13,:)=n1;  % x pixels
HdArray.Vals(14,:)=n2;  % y pixels

HdArray.Vals(61,:)=1;  % number of planes
HdArray.Vals(69,:)=33686018;  % machine type
end