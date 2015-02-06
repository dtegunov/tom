function m_out=tom_cov(m,flag)
%TOM_COV creates ...
%
%   m_out=tom_cov(m)
%
%PARAMETERS
%
%  INPUT
%   m                   ...
%  
%  OUTPUT
%   m_out       		...
%
%EXAMPLE
%   .. = tom_cov(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

if nargin <2
    flag='cov';
end;

error(nargchk(0, 2, nargin, 'struct'))

m_out=zeros(size(m,2),size(m,2),'single');

sz=sqrt(length(m(:,1)));

for i=1:size(m,2)
    for ii=(i):(size(m,2))
        if (strcmp(flag,'cov'))
           m_out(i,ii)=sum( (m(:,i)-(mean(m(:,i)))) .*(m(:,ii)-mean(m(:,ii))) )./ (size(m(:,i),1)-1);
        end;
        if (strcmp(flag,'corr1d'))  
            m_out(i,ii)=tom_ccc(m(:,i),m(:,ii),'norm'); 
        end;
        if (strcmp(flag,'corr2d'))
             m_out(i,ii)=tom_ccc(reshape(m(:,i),sz,sz),reshape(m(:,ii),sz,sz),'norm'); 
        end;
        
        if (strcmp(flag,'corr3d'))
             m_out(i,ii)=tom_ccc(reshape(m(:,i),sz,sz,sz),reshape(m(:,ii),sz,sz,sz),'norm'); 
        end;
        
        
    end;
    %i
end;

%fill up symmetric values
for i=2:size(m_out,1)
    for ii=1:(i-1)
        m_out(i,ii)=m_out(ii,i);
    end;
end;



%
% sz=80.*80;
%
% %n_of_lut=factorial(sz+1)./(2.*factorial(sz-1));
% n_of_lut = exp(gammaln(sz+2) - gammaln(3) - gammaln(sz-2));
% lut=zeros(n_of_lut,2,'single');
%
% zz=1;
% for i=1:sz
%     disp(i);
% %     for ii=(i):4
% %          lut(zz,1)=i;
% %          lut(zz,2)=ii;
% %         % disp(num2str([i ii]));
% %     end;
% %i
%     lut(zz:zz+sz-i,1)=i;
%     lut(zz:zz+sz-i,2)=i:sz;
%    % disp(num2str([i ii]));
%     zz=zz+sz-i+1;
% end;





% h=tom_reademheader(infile);
% sz_m=h.Header.Size;
%
% if nargin < 3
%     package=[1 sz_m(2)];
% end;
%
% if nargin < 4
%     buffer_size=5000;
% end;
%
% start=package(1);
% stop=package(2);


%if


% for i=1:sz_m(2)
%
%     out_buffer=zeros(sz_m(2),1);
%     m1=tom_emreadc(infile,'subregion',[1 i 1],[sz_m(1)-1 0 0]); m1=m1.Value;
%     m_tmp=m1;
%     buffer_count=buffer_size;
%     for ii=(i):sz_m(2)
%
%         if (buffer_count==buffer_size+1)
%             if (ii+buffer_size <= sz_m(2)+1)
%                 max_buf=buffer_size;
%             else
%                 max_buf=sz_m(2)-ii+1;
%             end;
%             m_tmp=tom_emreadc(infile,'subregion',[1 ii 1],[sz_m(1)-1 max_buf-1 0]); m_tmp=m_tmp.Value;
%             buffer_count=1;
%         end;
%
%         if (ii==i)
%             m2=m_tmp(:,1);
%         else
%             m2=m_tmp(:,buffer_count);
%         end;
%         out=sum( (m1-mean(m1)) .*(m2-mean(m2)) ) ./ (size(m1,1)-1);
%         out_buffer(ii)=out;
%         buffer_count=buffer_count+1;
%
%
%     end;
%
%     tom_emwritec(outfile,out_buffer','subregion',[i 1 1],[1 size(out_buffer,1) 1]);
%
%
% end;
%
%
%
% fill up symmetric values
% ii=1;
% for i=start+1:stop
%     buffer=tom_emreadc(outfile,'subregion',[ii i 1],[0 (sz_m(2)-i) 0]); buffer=buffer.Value;
%     tom_emwritec(outfile,buffer,'subregion',[i ii 1],[(sz_m(2)-i+1) 1 1]);
%     ii=ii+1;
% end;