function b=move(a,pix,cond);
%function b=move(a,pix,cond);
%processed image b,  input image a;
%pix=[x,y,z]; Move in x, y, z.
%cond: 0=Dirichlet, 1=reflecting, 2=periodic
%if cond = [] then Dirichlet Conditions

[s1,s2,s3]=size(a);

if nargin==2 cond=0; end

if s3==1;
  if cond==0 %Dirichlet boundary conditions
    if pix(1,1)<0
       b=[a(:,abs(pix(1,1))+1:s2),zeros(s1,abs(pix(1,1)))];
       a=b;
     elseif pix(1,1)>0
       b=[zeros(s1,abs(pix(1,1))),a(:,1:s2-abs(pix(1,1)))];
       a=b;
    end %if

    if pix(1,2)<0
       b=[a(abs(pix(1,2))+1:s1,:); zeros(abs(pix(1,2)),s2)];
       a=b;
     elseif pix(1,2)>0
       b=[zeros(abs(pix(1,2)),s2);a(1:s1-abs(pix(1,2)),:)];
       a=b;
    end %if

  elseif cond==1 %reflecting b.c.
    if pix(1,1)==-1
       b=[a(:,2:s2),a(:,s2)];
       a=b;
     elseif pix(1,1)==1
       b=[a(:,1),a(:,1:s2-1)];
       a=b;
    end %if

    if pix(1,2)==-1
       b=[a((2:s1),:);a(s1,:)];
       a=b;
     elseif pix(1,2)==1
       b=[a(1,:);a(1:s1-1,:)];
       a=b;
    end %if


  elseif cond==2 %periodic b.c.
    if pix(1,1)<0
       b=[a(:,abs(pix(1,1))+1:s2),a(:,1:abs(pix(1,1)))];
       a=b;
     elseif pix(1,1)>0
       b=[a(:,s2-abs(pix(1,1))+1:s2),a(:,1:s2-abs(pix(1,1)))];
       a=b;
    end %if

    if pix(1,2)<0
       b=[ a(abs(pix(1,2))+1:s1,:); a(1:abs(pix(1,2)),:) ];
       a=b;
     elseif pix(1,2)>0
       b=[ a(s1-abs(pix(1,2))+1:s1,:); a(1:s1-abs(pix(1,2)),:) ];
       a=b;
    end %if
  end %if

elseif s3>1;
  b=zeros(s1,s2,s3);
  if cond==0
    if pix(1,1)<0
       b(:,1:s1-abs(pix(1,1)),:)=a(:,abs(pix(1,1))+1:s2,:);
       a=b;
    elseif pix(1,1)>0
       b(:,abs(pix(1,1))+1:s2,:)=a(:,1:s2-abs(pix(1,1)),:);
       a=b;
    end %if

    if pix(1,2)<0
       b(1:s2-abs(pix(1,2)),:,:)=a(abs(pix(1,2))+1:s2,:,:);
       a=b;
    elseif pix(1,2)>0
       b(abs(pix(1,2))+1:s2,:,:)=a(1:s2-abs(pix(1,2)),:,:);      
       a=b;
    end %if

    if pix(1,3)<0
       b(:,:,1:s3-abs(pix(1,3)))=a(:,:,abs(pix(1,3))+1:s3);
       a=b;
    elseif pix(1,3)>0
       b(:,:,abs(pix(1,3))+1:s3)=a(:,:,1:s3-abs(pix(1,3)));
       a=b;
    end %if



  elseif cond==1
    if pix(1,1)<0
       b(:,1:s2-abs(pix(1,1)),:)=a(:,abs(pix(1,1))+1:s2,:);
       a=b;
    elseif pix(1,1)>0
       b(:,abs(pix(1,1))+1:s2,:)=a(:,1:s2-abs(pix(1,1)),:);
       a=b;
    end %if

    if pix(1,2)<0
       b=[a(abs(pix(1,2))+1:s1,:); zeros(abs(pix(1,2)),s2)];
       a=b;
    elseif pix(1,2)>0
       b=[zeros(abs(pix(1,2)),s2);a(1:s1-abs(pix(1,2)),:)];
       a=b;
    end %if
  end %if
end %if


