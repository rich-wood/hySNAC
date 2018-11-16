
%parameters

% dataset_name='SWE_EnvAccount' %this uses scb emission data, and not exiobase emissions data for sweden  - will use F and Fhh from excel files for sweden, and InputData\QEXIO for all other countries
dataset_name='chem' %chemicals - 'InputData\disagg_src\QEXIO_chem' - this uses the full stressor data provided by PRINCE consortium members in EXIOBASE classification already
% dataset_name='EXIOBASE_EnvAccount' %other exiobase stressors - will use InputData\QEXIO - this uses the full stressor data provided by PRINCE consortium members in EXIOBASE classification already

thisdir=pwd;

sweden_index=25




% prep exiobase data
% using directory: indecol\Projects\MRIOs\EXIOBASE3\EXIOBASE_3_3_20160218\Analysis\data
% this can be done in any other way, you just need to get the multipliers
% and the final demand out.
% cd d:\indecol\Projects\MRIOs\EXIOBASE3\EXIOBASE_3_3_20160218\Analysis\data


fl=1;flavor{2,fl}='pxp';
% for yr=2008:2014
for yr=2013
    yrstr=num2str(yr);
    cd('x:\indecol\Projects\MRIOs\EXIOBASE3\EXIOBASE_3_4\Analysis\data')
    
    [IO, meta] = load_mrio('EXIOBASE3', yrstr, 'EXIOBASE_3_4', 'pxp', 1);
    tr=mriotree(meta);
    
    if strcmp(dataset_name,'EXIOBASE_EnvAccount')
        load(['Qsrc_',yrstr,'_',flavor{2,fl},'.mat'])
%         load(['FsrcMult_',yrstr,'_',flavor{2,fl},'.mat'])
        load(['macro_',yrstr,'_',flavor{2,fl},'.mat'])
        
    elseif strcmp(dataset_name,'chem')
        
        load(['\\winfil.it.ntnu.no\EPT_eksperimentell\indecol\Projects\PRINCE\WP6\scripts\data\Qsrc',yrstr,'.mat'])
        load(['\\winfil.it.ntnu.no\EPT_eksperimentell\indecol\Projects\PRINCE\WP6\scripts\data\FsrcMult',yrstr,'.mat'])
        load(['\\winfil.it.ntnu.no\EPT_eksperimentell\indecol\Projects\PRINCE\WP6\scripts\data\macro',yrstr,'.mat'])
        
    end
    
    
    Qmult=Q10.mult.cnt; % these are the multipliersload('macro2011.mat')
    %     QEXIO.prod=Q10.prod.cnt(:,sweden_index);
    %     QEXIO.cons=Q10.cons.cnt(:,sweden_index);
    %     QEXIO.direct=Q10.direct.cnt(:,sweden_index);
    QEXIO.Yglobal=sum(macro.yall.sec.all,2); % final demand
    QEXIO.Yswe=macro.yall.sec.all(:,strcmp('SWE',meta.countries)); % final demand
    QEXIO.indic=indic;
    QEXIO.desc={'This is the multipliers and the final demand extracted from EIOXBASEv3.4, multipliers according to S*L and Y as the sum of all columns of the final demand block'}
    
    cd(thisdir)
    %%
    clear tmp*
    Cs=importdata(['..\..\Sweden Model\MetaData\MatchStressors_',dataset_name,'.xlsx'])
    for i=1:size(Cs.textdata,1)
        if ~isempty(Cs.textdata{i,3})
            fmatch=find(strcmp(Cs.textdata(i,3),QEXIO.indic.name));
            QEXIO.mult(:,i)=Qmult(fmatch(1),:);
            QEXIO.prod(:,i)=Q10.prod.cnt(fmatch(1),sweden_index);
            QEXIO.prod_ind(:,i)=Q10.prod_ind.sec(fmatch(1),:,sweden_index);
            QEXIO.cons(:,i)=Q10.cons.cnt(fmatch(1),sweden_index);
            QEXIO.direct(:,i)=Q10.direct.cnt(fmatch(1),sweden_index);
            QEXIO.mult_src_cnt(:,:,i)=sparse(tr.collapseZdim(diag(sum(Q10.int.cnt(fmatch(1),:),1)),1))*IO.L;
            QEXIO.mult_src_sec(:,:,i)=tr.sectorSum(sparse(diag(sum(Q10.int.cnt(fmatch(1),:),1)))*IO.L,1);
        else
            for j=1:Cs.data(i,2)
                fmatch=find(strncmp(Cs.textdata{i,6+j},meta.labsF,length(Cs.textdata{i,6+j})));
                if j==1
                    tmp_src=sparse(diag(sum(IO.S(fmatch,:),1)))*IO.L;
                    tmp=sum(Fm.mult.cnt(fmatch,:),1);
                    tmp1=sum(F.prod.cnt(fmatch,sweden_index),1);
                    tmp11=sum(F.prod_ind.sec(fmatch,sweden_index),:,1);
                    tmp2=sum(F.cons.cnt(fmatch,sweden_index),1);
                    tmp3=sum(F.direct.cnt(fmatch,sweden_index),1);
                else
                    tmp_src=tmp_src+sparse(diag(sum(IO.S(fmatch,:),1)))*IO.L;
                    tmp=tmp+sum(Fm.mult.cnt(fmatch,:),1);
                    tmp1=tmp1+sum(F.prod.cnt(fmatch,sweden_index),1);
                    tmp11=tmp1+sum(F.prod_ind.sec(fmatch,:,sweden_index),1);
                    tmp2=tmp2+sum(F.cons.cnt(fmatch,sweden_index),1);
                    tmp3=tmp3+sum(F.direct.cnt(fmatch,sweden_index),1);
                end
            end
            QEXIO.mult(:,i)=tmp;
            QEXIO.prod(:,i)=tmp1;
            QEXIO.prod_ind(:,i)=tmp1;
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
        QEXIO.prod_ind(:,i)=QEXIO.prod_ind(:,i)/Cs.data(i,1);
        QEXIO.cons(:,i)=QEXIO.cons(:,i)/Cs.data(i,1);
        QEXIO.direct(:,i)=QEXIO.direct(:,i)/Cs.data(i,1);
    end
    save(['QEXIO_',dataset_name,yrstr],'QEXIO','meta','indic')
end