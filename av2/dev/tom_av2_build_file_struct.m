function tom_av2_build_file_struct(flag)


if nargin==0
    flag='';
end;

%rock the struct

mkdir align;
    mkdir align/high/;
        mkdir align/high/whole_particle/;
        mkdir align/high/core/;
    mkdir align/low/;
        mkdir align/low/whole_particle/;
        mkdir align/low/core/;

mkdir pick;
    mkdir pick/high/;
    mkdir pick/low/;
        
mkdir org;
    
mkdir sort;
    mkdir sort/high/;
    mkdir sort/low/; 

mkdir rec1;
    mkdir rec1/high/;
        mkdir rec1/high/step1;
        mkdir rec1/high/step2;
        mkdir rec1/high/step3;
        mkdir rec1/high/step4;
        mkdir rec1/high/step5;
    mkdir rec1/low/; 
        mkdir rec1/low/step1;
        mkdir rec1/low/step2;
        mkdir rec1/low/step3;
        mkdir rec1/low/step4;
        mkdir rec1/low/step5;
        
    
mkdir obsolete;
