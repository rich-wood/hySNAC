% SWEDEN_FIRST_IO_Model.m - Main model file that takes pre-prepared IO
% multipliers from a MRIO and links it to a country-level IO table.

% Manuscript Title: Hybrid SNAC for calculation of environmental footprints � using life-cycle approaches via input-output multipliers on traded goods
% Contact: Richard Wood
% richard.wood@ntnu.no


% Master script: Current

% Dependencies: Pre-calculated data in InputData. These can be called by
% this scripts by setting prepdata_for_sweden to 1.

% Addtional comments:



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define parameters here:

% year_of_analysis=2008:2014
year_of_analysis=2014

dataset_name='SWE' %this uses scb emission data, and not exiobase emissions data for sweden  - will use F and Fhh from excel files for sweden, and InputData\QEXIO for all other countries
% dataset_name='chem' %chemicals - 'InputData\disagg_src\QEXIO_chem' - this uses the full stressor data provided by PRINCE consortium members in EXIOBASE classification already
% dataset_name='EXIOBASE_EnvAccount' %other exiobase stressors - will use InputData\QEXIO - this uses the full stressor data provided by PRINCE consortium members in EXIOBASE classification already

%
update_results_template=0
prepdata_for_sweden=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Other parameters that should be automatically assigned:

% setting up input files:
exiobase_input_data_path=['InputData\QEXIO_',dataset_name];

stressor_raw.name=readtable(['MetaData\MatchStressors_',dataset_name,'.xlsx'],'Sheet','Names');
stressor_raw.Concordance=readtable(['MetaData\MatchStressors_',dataset_name,'.xlsx'],'Sheet','Concordance');
indic_hy.name=stressor_raw.name.LongName;
indic_hy.sh_name=stressor_raw.name.ShortName;
indic_hy.unit=stressor_raw.name.Unit;
indic_hy.src=stressor_raw.name.Source;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%    %prepdata
if prepdata_for_sweden==1
    cd Inputdata
    for yr=year_of_analysis
        prep_swe_data(yr)
%         prep_exio_data_list(yr,indic_hy,stressor_raw)
    end
    cd ..
        
end

%%
if update_results_template
    for yr=year_of_analysis
        yrstr=num2str(yr);
        resultfile=['Results\SWE_SRIO_',dataset_name,'_results_',yrstr,'.xlsx']
        copyfile('Results\template.xlsx',resultfile)
    end
end


%%

    %aggregation matrix to prod
    Conc_prod_data=importdata('MetaData\Agg_ExIOB3_MirSNI07_MB_agg65_59.xlsx');
    C_ind=Conc_prod_data.data.Key59(2:60,4:end);
    C_prod_to_ind=importdata('MetaData\EXIOBASE20p_EXIOBASE20i.txt');
    C_prod=C_ind*C_prod_to_ind.data';

%%
for yr=year_of_analysis
    yrstr=num2str(yr);
    
    resultfile=['Results\SWE_SRIO_',dataset_name,'_results_',yrstr,'.xlsx']
    
    
    % Load Sweden data:
    load(['InputData\SWE',yrstr],'SWE')
    % Load exiobase data:
    load([exiobase_input_data_path,yrstr],'QEXIO','meta')
    
    %% BEGIN CLASSIFICATION MATCHING
    
    %country names:
    if exist('MetaData\CountryMapping.mat')
        load('MetaData\CountryMapping')
    else
        cnt_names=SWE.textdata.(['x',yrstr])(219,3:201+3);
        C_cnt=indCountry(cnt_names);
        save('MetaData\CountryMapping','C_cnt')
    end
    %assume RoW is region 45
    C_cnt.regi(end)=45;
    
    C_prod_multi_region=kron(eye(meta.NCOUNTRIES),C_prod'); %G^p in paper
    
    % Calculate multipliers
    cd Inputdata
    [F,Fhh]=load_compile_stressors(indic_hy,yr,QEXIO,SWE,C_prod);
    cd ..
    
    Sdom_all=F/diag(SWE.x); %environmental extension intensities
    
    %these are precalculated multipliers from EXIOBASE3.4 for a number of environmental impacts (named in the indic variable)
%          Multipliers for all extensions from exiobase
    IOmult=QEXIO.mult';
    IOmult_src=QEXIO.mult_src;
    
    n_indic=length(indic_hy.name);
    
    % now calculate the total impact by product of the aggregated products:
    % (multiper*demand)*aggregation matrix
    % SLy ?G^p - (SLy ?, where y is global final demand by product for EXIOBASE ), aggregating to Swedish classification by post-multiplying by G^p
    tmp_Total_impact=(IOmult*diag(QEXIO.Yglobal))*C_prod_multi_region;
    % do same thing for total final demand from EXIOBASE "aggregation of
    % global final demand (yG^p)  "
    Total_Y=(QEXIO.Yglobal)'*C_prod_multi_region;
    
    
    %Repeat but with source detailed:
    tmp_y=(sparse(diag(QEXIO.Yglobal)))*sparse(C_prod_multi_region);
    Total_impact_src=zeros(meta.NCOUNTRIES,size(tmp_y,2),n_indic);
    for i=1:n_indic
        Total_impact_src(:,:,i)=IOmult_src(:,:,i)*tmp_y;
    end
    
    
    %Exchange rate:
    exchrate=9.09;
    Total_Y_SEK=Total_Y*exchrate;
    % reverse engineer the aggregate multipliers by dividing the aggregated
    % impact by product group and region by the final demand by product group
    % and region: 
%   Q^g=SLy ?G^p*  ?((yG^p ) ?)?^(-1)
%   Q^g has thus a column dimension of 49 (EXIOBASE countries) x 59 (Swedish products)
    EXIO_Agg_multipliers_SEK=tmp_Total_impact./repmat(eps+Total_Y_SEK,size(tmp_Total_impact,1),1);
    EXIO_Agg_multipliers_SEK_src=Total_impact_src./(eps+(repmat(Total_Y_SEK,[meta.NCOUNTRIES,1,n_indic])));
    
    
    %To do if not 2008-2014:
    % price indices
    
    clear footprint_by_src_byprod footprint_by_src
    % Connect to Excel
    Excel = actxserver('excel.application');
    % Get Workbook object
    Excel.Workbooks.Open([pwd,'\',resultfile]);
    
    %% loop over indicators
    for i_indic=1:n_indic
        
        
        
        EXIO_Agg_multipliers_cnt_by_prod=reshape(EXIO_Agg_multipliers_SEK(i_indic,:),59,meta.NCOUNTRIES)';
        
        
        EXIO_Agg_multipliers_cnt_by_prod_by_src_cnt.val=reshape(EXIO_Agg_multipliers_SEK_src(:,:,i_indic),meta.NCOUNTRIES,59,meta.NCOUNTRIES);
        EXIO_Agg_multipliers_cnt_by_prod_by_src_cnt.dim1={'Country of origin EXIOBASE'};
        EXIO_Agg_multipliers_cnt_by_prod_by_src_cnt.dim2={'Product'};
        EXIO_Agg_multipliers_cnt_by_prod_by_src_cnt.dim3={'Country of import EXIOBASE'};
        
        %disaggregate to the SWEDEN country classification (circa 202
        %countries). Asssumes equal multipliers of an aggregate EXIO region
        %to a country in that region.
%         in the paper: G^c is a binary matrix  
%       Here instead of using G^c, we use indexing (C_cnt.reg), as it is soo much
%       easier to code
        EXIO_Agg_multipliers_disagg_cnt_by_prod=EXIO_Agg_multipliers_cnt_by_prod(C_cnt.regi,:);
        
        
        EXIO_Agg_multipliers_disagg_cnt_by_prod_by_src_cnt.val=EXIO_Agg_multipliers_cnt_by_prod_by_src_cnt.val(:,:,C_cnt.regi);
        EXIO_Agg_multipliers_disagg_cnt_by_prod_by_src_cnt.dim1={'Country of origin EXIOBASE'};
        EXIO_Agg_multipliers_disagg_cnt_by_prod_by_src_cnt.dim2={'Product'};
        EXIO_Agg_multipliers_disagg_cnt_by_prod_by_src_cnt.dim3={'Country of import SWE'};
        
        % Now we call our final regionally disaggregated impact per unit of sweden imports "Qfin"
        % These are the multipliers x by the trade shares (which add up to
        % 1). Hence Qfin are not multipliers, per-se, but multipliers
        % disaggregated by country of final import
        %
        %if we aggregate country of origin (of env impact):
        %         Q_fin=sum(Agg_multipliers_disagg_cnt_by_prod.*trade_shares',1);
        % or we can do this with the country of origin resolution
        Q_fin=EXIO_Agg_multipliers_disagg_cnt_by_prod.*SWE.trade_shares';
        
        for i_cnt=1:meta.NCOUNTRIES
            Q_fin_src_cnt.val(i_cnt,:,:)=squeeze(EXIO_Agg_multipliers_disagg_cnt_by_prod_by_src_cnt.val(i_cnt,:,:)).*SWE.trade_shares;
            Q_fin_src_cnt.dim1={'Country of origin EXIOBASE'};
            Q_fin_src_cnt.dim2={'Product'};
            Q_fin_src_cnt.dim3={'Country of import SWE'};
        end
        
        % For Sweden, extract relevant indicator:
        Sdom=Sdom_all(i_indic,:);
        
        %Now we can calculate the full model as a result of the domestic component
        %+ the impacts embodied in intermediate imports + the impacts embodied in
        %final imports:
        SwedenFP_dom_imp_incl_exports.meta={'domestic component + the impacts embodied in intermediate imports + the impacts embodied in final imports:'};
        
        for j=1:8 % loop over final demand components (incl exports=column8)
            SwedenFP_dom_imp_incl_exports.val(:,j)=Sdom*SWE.Ldom*diag(SWE.Yd(:,j))+sum(Q_fin*SWE.Aimp*SWE.Ldom,1)*diag(SWE.Yd(:,j))+sum(Q_fin,1)*diag(SWE.Ym(:,j));
        end
        
        %Now exclude exports from swedish demand:
        SWE.Yd_dom=sum(SWE.Yd(:,1:7),2);
        SWE.Ym_dom=sum(SWE.Ym(:,1:7),2);
        %Sweden footprint
        SWE_cons(i_indic)=Sdom*SWE.Ldom*SWE.Yd_dom...
            +sum(Q_fin*SWE.Aimp*SWE.Ldom*SWE.Yd_dom,1)...
            +sum(Q_fin*SWE.Ym_dom,1)...
            +sum(Fhh(i_indic,:));
        % Sweden footprint disagg by product:
        SWE_cons_by_prod(:,i_indic)=(Sdom*SWE.Ldom*diag(SWE.Yd_dom)...
            +sum(Q_fin*SWE.Aimp*SWE.Ldom,1)*diag(SWE.Yd_dom)...
            +sum(Q_fin,1)*diag(SWE.Ym_dom))';
        
        % Exports is column 8, calculate them too:
        SWE_export_by_prod(:,i_indic)=(Sdom*SWE.Ldom*diag(SWE.Yd(:,8))...
            +sum(Q_fin*SWE.Aimp*SWE.Ldom,1)*diag(SWE.Yd(:,8))...
            +sum(Q_fin,1)*diag(SWE.Ym(:,8)))';
        
        
        SWE_production(:,i_indic)=(Sdom*diag(SWE.x));
        
        %for interest, not really used:
        SwedenFP_dom_imp_excl_exports_by_prod_dom_only(:,i_indic)=(Sdom*SWE.Ldom*diag(SWE.Yd_dom))';
        
        %for interest, not really used:
        compare_mult=[(Sdom*SWE.Ldom)',Q_fin'];
        
        %Now also calculate emissions embodied in all imports (by product)
        imported_emissions=(Q_fin*diag(SWE.M))';
        imported_emissions_by_indic(:,i_indic)=sum(imported_emissions,1)';
        imported_emissions_by_indic_by_prod(:,i_indic)=sum(imported_emissions,2);
        
        %
        %         imported emissions by source using Q_fin_src_cnt
        for j_cnt=1:49
            imported_emissions_by_src(j_cnt,i_indic)=sum(sum((squeeze(Q_fin_src_cnt.val(j_cnt,:,:))'*diag(SWE.M))',1),2);
        end
        
        %         need to implement emissions to domestic final demand by source using Q_fin_src_cnt
        for j_cnt=1:49
            footprint_by_src_byprod(j_cnt,:,i_indic)=...            
            +sum(squeeze(Q_fin_src_cnt.val(j_cnt,:,:))'*SWE.Aimp*SWE.Ldom,1)*diag(SWE.Yd_dom)...
            +sum(squeeze(Q_fin_src_cnt.val(j_cnt,:,:))',1)*diag(SWE.Ym_dom);
        end
        
        footprint_by_src_byprod(25,:,i_indic)=footprint_by_src_byprod(25,:,i_indic)+...        
        Sdom*SWE.Ldom*diag(SWE.Yd_dom);
        
        
        %compare EXIOBASE to Sweden Model:
        compare_prod_account(i_indic,:)=...
            [QEXIO.prod(i_indic),sum(F(i_indic,:),2)+sum(Fhh(i_indic,:),2),QEXIO.prod(i_indic)/(sum(F(i_indic,:),2)+sum(Fhh(i_indic,:),2))];
        compare_cons_account(i_indic,:)=...
            [QEXIO.cons(i_indic),SWE_cons(i_indic),QEXIO.cons(i_indic)/SWE_cons(i_indic)];
        compare_direct(i_indic,:)=...
            [QEXIO.F_hh(i_indic),sum(Fhh(i_indic,:),2),QEXIO.F_hh(i_indic)/sum(Fhh(i_indic,:),2)];
        
        
        short_indic_name=indic_hy.sh_name{i_indic};
        
          %%  Begin Results writing....somewhat messy!
        if n_indic<20
        
        % write all results by country
        xlswrite1(resultfile,[SwedenFP_dom_imp_incl_exports.val,sum(imported_emissions,2)],short_indic_name,'c3')
%         xlswrite1(resultfile,sum([SwedenFP_dom_imp_incl_exports.val,sum(imported_emissions,2)],1),short_indic_name,'c2')
        xlswrite1(resultfile,[SWE.Fd_names,'Imports'],short_indic_name,'c1')
        xlswrite1(resultfile,[SWE.prod_names;'Imports'],short_indic_name,'a3')
%         xlswrite1(resultfile,{'Total'},short_indic_name,'a2')
        xlswrite1(resultfile,{'Household'},short_indic_name,'a63')
        xlswrite1(resultfile,Fhh(i_indic,:),short_indic_name,'c63')
        
        %if you want to write country of final production/export results,
        %uncomment this:
%         xlswrite1(resultfile,[imported_emissions],[short_indic_name,'_imp'],'c4')
%         xlswrite1(resultfile,SWE.cnt_names,[short_indic_name,'_imp'],'c1')
%         xlswrite1(resultfile,SWE.prod_names,[short_indic_name,'_imp'],'a4')
%         
        
        %if you want to write country of origin by product results,
        xlswrite1(resultfile,footprint_by_src_byprod(:,:,i_indic)',[short_indic_name,'_origin'],'c4')
        xlswrite1(resultfile,meta.countrynames',[short_indic_name,'_origin'],'c1')
        xlswrite1(resultfile,SWE.prod_names,[short_indic_name,'_origin'],'a4')
        xlswrite1(resultfile,{'Household'},[short_indic_name,'_origin'],'a63')
        originFhh=zeros(1,49);
        originFhh(25)=Fhh(i_indic,:);
        xlswrite1(resultfile,originFhh,[short_indic_name,'_origin'],'c63')
%         

        end
        
    end  % looping over indicators
    
    
    footprint_by_src=squeeze(sum(footprint_by_src_byprod,2));
    footprint_by_src(25,:)=footprint_by_src(25,:)+Fhh';
    
    
    
    %% Results writing....somewhat messy!
    
    xlswrite1(resultfile,compare_prod_account,'EXIO3vsSWE','c3')
    xlswrite1(resultfile,indic_hy.name,'EXIO3vsSWE','a3')
    xlswrite1(resultfile,indic_hy.unit,'EXIO3vsSWE','b3')
    xlswrite1(resultfile,compare_cons_account,'EXIO3vsSWE','c23')
    xlswrite1(resultfile,indic_hy.name,'EXIO3vsSWE','a23')
    xlswrite1(resultfile,indic_hy.unit,'EXIO3vsSWE','b23')
    xlswrite1(resultfile,compare_direct,'EXIO3vsSWE','c43')
    xlswrite1(resultfile,indic_hy.name,'EXIO3vsSWE','a43')
    xlswrite1(resultfile,indic_hy.unit,'EXIO3vsSWE','b43')
    xlswrite1(resultfile,SWE_production,'Production','c3')
    xlswrite1(resultfile,sum(SWE_production),'Production','c62')
    xlswrite1(resultfile,Fhh','Production','c63')
    xlswrite1(resultfile,sum(SWE_production)+Fhh','Production','c64')
    xlswrite1(resultfile,indic_hy.name','Production','c1')
    xlswrite1(resultfile,indic_hy.unit','Production','c2')
    xlswrite1(resultfile,SWE_cons_by_prod,'Cons_Product','c3')
    xlswrite1(resultfile,sum(SWE_cons_by_prod),'Cons_Product','c62')
    xlswrite1(resultfile,Fhh','Cons_Product','c63')
    xlswrite1(resultfile,sum(SWE_cons_by_prod)+Fhh','Cons_Product','c64')
    xlswrite1(resultfile,indic_hy.name','Cons_Product','c1')
    xlswrite1(resultfile,indic_hy.unit','Cons_Product','c2')
    xlswrite1(resultfile,imported_emissions_by_indic_by_prod,'Imports','c3')
    xlswrite1(resultfile,indic_hy.name','Imports','c1')
    xlswrite1(resultfile,indic_hy.unit','Imports','c2')
    xlswrite1(resultfile,SWE_export_by_prod,'Exports','c3')
    xlswrite1(resultfile,indic_hy.name','Exports','c1')
    xlswrite1(resultfile,indic_hy.unit','Exports','c2')
    
    
    xlswrite1(resultfile,SWE_cons_by_prod./(SWE.Yd_dom+SWE.Ym_dom),'Cons_Product_intensity','c3')
    xlswrite1(resultfile,indic_hy.name','Cons_Product_intensity','c1')
    xlswrite1(resultfile,indic_hy.unit','Cons_Product_intensity','c2')
    xlswrite1(resultfile,imported_emissions_by_indic_by_prod./repmat(SWE.M,1,n_indic),'Imports_intensity','c3')
    xlswrite1(resultfile,indic_hy.name','Imports_intensity','c1')
    xlswrite1(resultfile,indic_hy.unit','Imports_intensity','c2')
    xlswrite1(resultfile,SWE_export_by_prod./repmat((SWE.Yd(:,8)+SWE.Ym(:,8)),1,n_indic),'Exports_intensity','c3')
    xlswrite1(resultfile,indic_hy.name','Exports_intensity','c1')
    xlswrite1(resultfile,indic_hy.unit','Exports_intensity','c2')
    
    
    xlswrite1(resultfile,imported_emissions_by_indic,'Country_Export','c4')
    xlswrite1(resultfile,SWE.cnt_names','Country_Export','a4')
    xlswrite1(resultfile,indic_hy.name','Country_Export','c1')
    xlswrite1(resultfile,indic_hy.unit','Country_Export','c2')
    
    
    
    xlswrite1(resultfile,[footprint_by_src],['Footprint_Origin'],'c4')
    xlswrite1(resultfile,sum([footprint_by_src]),['Footprint_Origin'],'c53')
    xlswrite1(resultfile,meta.countrynames,['Footprint_Origin'],'a4')
    xlswrite1(resultfile,indic_hy.name','Footprint_Origin','c1')
    xlswrite1(resultfile,indic_hy.unit','Footprint_Origin','c2')
    
    xlswrite1(resultfile,[imported_emissions_by_src],['Total_Imp_by_Country'],'c4')
    xlswrite1(resultfile,sum([imported_emissions_by_src]),['Total_Imp_by_Country'],'c53')
    xlswrite1(resultfile,meta.countrynames,['Total_Imp_by_Country'],'a4')
    xlswrite1(resultfile,indic_hy.name','Total_Imp_by_Country','c1')
    xlswrite1(resultfile,indic_hy.unit','Total_Imp_by_Country','c2')
    
    
    
    
    Excel.ActiveWorkbook.Save; Excel.ActiveWorkbook.Close; Excel.Quit; clear Excel
    disp('finished')
    
end
