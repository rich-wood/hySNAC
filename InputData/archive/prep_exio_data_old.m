% prep exiobase data
% using directory: indecol\Projects\MRIOs\EXIOBASE3\EXIOBASE_3_3_20160218\Analysis\data
% this can be done in any other way, you just need to get the multipliers
% and the final demand out.
% cd d:\indecol\Projects\MRIOs\EXIOBASE3\EXIOBASE_3_3_20160218\Analysis\data
fl=1;flavor{2,fl}='pxp';
for yr=2008:2014
    yrstr=num2str(yr);
    cd('x:\indecol\Projects\MRIOs\EXIOBASE3\EXIOBASE_3_4\Analysis\data')
    
    [IO, meta] = load_mrio('EXIOBASE3', yr, 'EXIOBASE_3_4', flavor{2,fl}, 1);
    tr=mriotree(meta);
    
    load(['Qsrc_',flavor{2,fl},yrstr,'.mat'])
    load(['FsrcMult_',flavor{2,fl},yrstr,'.mat'])
    load(['macro_',flavor{2,fl},yrstr,'.mat'])
    Qmult=Q10.mult.cnt; % these are the multipliersload('macro2011.mat')
%     QEXIO.prod=Q10.prod.cnt(:,25);
%     QEXIO.cons=Q10.cons.cnt(:,25);
%     QEXIO.direct=Q10.direct.cnt(:,25);
    QEXIO.Yglobal=sum(macro.yall.sec.all,2); % final demand
    QEXIO.Yswe=macro.yall.sec.all(:,strcmp('SWE',meta.countries)); % final demand
    QEXIO.indic=indic;
    QEXIO.desc={'This is the multipliers and the final demand extracted from EIOXBASEv3.3, multipliers according to S*L and Y as the sum of all columns of the final demand block'}
    cd('C:\Users\richardw\Box Sync\Projects\PRINCE_Swedish_EPA\WP6\Sweden Model\InputData\disagg_src')
    %%
    clear tmp*
    Cs=importdata('C:\Users\richardw\Box Sync\Projects\PRINCE_Swedish_EPA\WP6\Sweden Model\MetaData\MatchStressors.xlsx')
    for i=1:14
        if ~isempty(Cs.textdata{i,3})
            fmatch=find(strcmp(Cs.textdata(i,3),QEXIO.indic.name));
            QEXIO.mult(:,i)=Qmult(fmatch(1),:);
            QEXIO.prod(:,i)=Q10.prod.cnt(fmatch(1),25);
            QEXIO.cons(:,i)=Q10.cons.cnt(fmatch(1),25);
            QEXIO.direct(:,i)=Q10.direct.cnt(fmatch(1),25);
            QEXIO.mult_src_cnt(:,:,i)=sparse(tr.collapseZdim(diag(sum(Q10.int.cnt(fmatch(1),:),1)),1))*IO.L;
            QEXIO.mult_src_sec(:,:,i)=tr.sectorSum(sparse(diag(sum(Q10.int.cnt(fmatch(1),:),1)))*IO.L,1);
        else
            for j=1:Cs.data(i,2)
                fmatch=find(strncmp(Cs.textdata{i,6+j},meta.labsF,length(Cs.textdata{i,6+j})));
                if j==1
                    tmp_src=sparse(diag(sum(IO.S(fmatch,:),1)))*IO.L;
                    tmp=sum(Fm.mult.cnt(fmatch,:),1);
                    tmp1=sum(F.prod.cnt(fmatch,25),1);
                    tmp2=sum(F.cons.cnt(fmatch,25),1);
                    tmp3=sum(F.direct.cnt(fmatch,25),1);
                else    
                    tmp_src=tmp_src+sparse(diag(sum(IO.S(fmatch,:),1)))*IO.L;
                    tmp=tmp+sum(Fm.mult.cnt(fmatch,:),1);
                    tmp1=tmp1+sum(F.prod.cnt(fmatch,25),1);
                    tmp2=tmp2+sum(F.cons.cnt(fmatch,25),1);
                    tmp3=tmp3+sum(F.direct.cnt(fmatch,25),1);
                end
            end
            QEXIO.mult(:,i)=tmp;
            QEXIO.prod(:,i)=tmp1;
            QEXIO.cons(:,i)=tmp2;
            QEXIO.direct(:,i)=tmp3;
            QEXIO.mult_src_sec(:,:,i)=tr.sectorSum(tmp_src,1);
            QEXIO.mult_src_cnt(:,:,i)=tr.collapseZdim(tmp_src,1);
            
        end
        
        %unit conversion:
        QEXIO.mult(:,i)=QEXIO.mult(:,i)/Cs.data(i,1);
        QEXIO.mult_src_cnt(:,i)=QEXIO.mult_src_cnt(:,i)/Cs.data(i,1);
        QEXIO.mult_src_sec(:,i)=QEXIO.mult_src_sec(:,i)/Cs.data(i,1);
        QEXIO.prod(:,i)=QEXIO.prod(:,i)/Cs.data(i,1);
        QEXIO.cons(:,i)=QEXIO.cons(:,i)/Cs.data(i,1);
        QEXIO.direct(:,i)=QEXIO.direct(:,i)/Cs.data(i,1);
    end
    save(['QEXIO',yrstr],'QEXIO')
end