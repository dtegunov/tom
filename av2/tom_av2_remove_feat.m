function tom_av2_remove_feat(in_sel,out_sel,tmpl,perm_mask,peak_rad,thresh,max_num)
%TOM_AV2_REMOVE_FEAT premutes background pixels
%   
%
%  tom_av2_remove_feat(doc_filename,f_mask3d)
%
%  TOM_AV2_REMOVE_FEAT premutes background pixels
%  
%
%PARAMETERS
%
%  INPUT
%   doc_filename      *.doc filename 
%   f_mask3d          3d mask filename
%   
%
%  OUTPUT
%
%EXAMPLE
%     
%  tom_av2_permute_bg('all_cut.doc','mask3d.em');
%
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_em_classify3d
%
%   created by FB 08/09/09
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

warning off;

tt=importdata('all64_org.sel');
img=tom_spiderread(tt.textdata{1});
sz=size(img.Value);
mid=floor(sz./2)+1;
figure; 
for i=1:length(tt.textdata)
    img=tom_spiderread(tt.textdata{i});
    ime=mean2(img.Value);
    istd=std2(img.Value);
    
    cc=tom_os3_corr(img.Value,tmpl,'norm');
    tom_imagesc(img.Value);
    cc_fun=cc;
    for ii=1:max_num
        [pos(ii,:)  val(ii) cc_fun]=tom_peak(cc_fun,peak_rad);
    end
    
    idx=find(val>0.56);
    
    img_clean=img.Value;
    for ii=1:length(idx)
        sh=pos(idx(ii),:)-mid;
        m_tmp=tom_move(perm_mask,sh);
        rnd_str=tom_norm(rand(length(find(m_tmp>0)),1),'mean0+1std');
        img_clean(find(m_tmp>0))=((rnd_str.*istd)+ime);
        %img_perm=tom_permute_bg(img.Value,m_tmp==0);
    end;
    subplot(1,3,1); tom_imagesc(img.Value); hold on; plot(pos(idx(:),1),pos(idx(:),2),'ro'); hold off;
    subplot(1,3,2); plot(val,'r-'); 
    subplot(1,3,3); tom_imagesc(img_clean); 
    
    
    disp(' ');
end;
