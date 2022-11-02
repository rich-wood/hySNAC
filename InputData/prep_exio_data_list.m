%OPTIONAL Function% function prep_exio_data(yr,indic_hy,stressor_raw)

%The purpose of this function is to prep MRIO data for use in the Sweden
%hybrid model.
%It produces all required intensities and multipliers from the MRIO model,
%that can then be used in a single country model (i.e. it handles all the
%computationally heavy MRIO calculations).
%It produces a number of "indicators" which may directly correspond to
%substances, or be a characterisation of a number of substances.
% Key information is
%       1. MRIO data (default EXIOBASE3)
%       2. Stressors/environmental intensities and matching meta data (if
%       not default from EXIOBASE, must be specified as an additional
%       "dataset"
%       3. A meta correspondence file that matches the meta data associated
%       with the environmental intensities with the indicators.

% Key output is a structure QEXIO with the data, meta, with updated
% metainformation, a country index for sweden corresponding to the data and
% a structure "matched_stressors" which shows the final relationship
% between stressors and indicators.

%% parameters

%You must choose a dataset here that corresponds to fielnames.
% dataset_name='chem' %chemicals -  - this uses the full stressor data provided by PRINCE consortium members in EXIOBASE classification already - see line 68 for updating file names
dataset_name='SWE' % Sweden account
ntnu_user=0; % choose if you are on the NTNU network, otherwise you will load files from:  https://ntnu.box.com/s/75v6quxvuxlpw7fu4roa7lnrfc09jotr
save_data_for_others=0; %if you are an NTNU user and want to update the saved files, change flag to 1.
load_io_data=1 %flag to cut run time if data pre-loaded - be careful on which year you have data pre-loaded for!!!



if strcmp(dataset_name,'SWE')
    year_of_analysis=2008:2014 
else
    year_of_analysis=2013 %only 2013 for chemical paper
end


if length(year_of_analysis)>1
    load_io_data=1
end
if ~exist('meta')
    load_io_data=1
end

thisdir=pwd;


%% load meta information about relationship between indicators and substances:
stressor_raw.name=readtable(['..\MetaData\MatchStressors_',dataset_name,'.xlsx'],'Sheet','Names');
stressor_raw.Concordance=readtable(['..\MetaData\MatchStressors_',dataset_name,'.xlsx'],'Sheet','Concordance');


indic_hy.name=stressor_raw.name.LongName;
indic_hy.sh_name=stressor_raw.name.ShortName;
indic_hy.unit=stressor_raw.name.Unit;
indic_hy.src=stressor_raw.name.Source;



%% prep exiobase data
% use load function and extract intensities and labels.

for yr=year_of_analysis
    
    yrstr=num2str(yr);
    
    
    if ntnu_user==1
        
        addpath('\\winfil.it.ntnu.no\EPT_eksperimentell\indecol\Software\GeneralMatlabUtilities\')
        
        if load_io_data==1
            fl='pxp'
            [IO, meta] = load_mrio('EXIOBASE3', yrstr, 'EXIOBASE_3_41', fl, 1);
            load(['\\winfil.it.ntnu.no\EPT_eksperimentell\indecol\Projects\MRIOs\EXIOBASE3\EXIOBASE_3_41\Analysis\data\macro_',yrstr,'_',fl,'.mat'],'macro')
        end
        
        if strcmp(dataset_name,'SWE')
            % No need to load, just using default EXIOBASE data
            stressors=IO.S;
            stressors_hh=IO.F_hh;
            
        elseif strcmp(dataset_name,'chem') % CHECK HERE IT PROBABLY NEEDS UPDATING %
            load(['\\winfil.it.ntnu.no\EPT_eksperimentell\indecol\Projects\PRINCE\WP6\processedData\extensions_PRINCE_',yrstr,'_chemPaper.mat'])
            if strcmp(fl,'ixi')
                stressors=[S_PRINCE_i;IO.S];
            else
                stressors=[S_PRINCE_p;IO.S];
            end
            stressors_hh=[zeros(size(S_PRINCE_p,1),size(IO.F_hh,2));IO.F_hh];
            meta.labsF=[labsF_PRINCE;IO.meta.labsF(:,1:3)];
            
            %         load('chem_tox_water_air.mat')
            %         S_tox_chem_tmp=tox.val/(sparse(diag(IO.x+double(IO.x<1e-1))));
            %         stressors=[S_PRINCE_i(1:2,:);S_tox_chem_tmp;S_PRINCE_i(137:140,:)];
            %         stressors_hh=zeros(size(stressors,1),1);
            %         meta.labsF=[labsF_PRINCE(1:2,:);tox.names;labsF_PRINCE(137:140,:)];
            
        end
        
        %
        %%%% if you want to save the input data for others to use:
        if save_data_for_others==1
            IO=rmfield(IO,'A');
            IO=rmfield(IO,'Z');
            cd 'C:\Users\richardw\Box Sync\Projects\Prince_Data'
            save(['prince_data_',dataset_name,'_',fl,'_',yrstr],'IO','meta','macro','stressors','stressors_hh')
        end
    else
        fl='pxp'
        %NB - you must have downloaded a local copy of the files from https://ntnu.box.com/s/75v6quxvuxlpw7fu4roa7lnrfc09jotr
%         load(['prince_data_',dataset_name,'_',fl,'_',yrstr],'IO','meta','macro','stressors','stressors_hh')
            load(['..\..\..\..\..\Prince_Data\prince_data_',dataset_name,'_',fl,'_',yrstr],'IO','meta','macro','stressors','stressors_hh')
            addpath('..\..\..\..\..\Prince_Data\')
    end
    
    %Final prep of data:
    
%     tr=mriotree(meta);
    sweden_index=25;%tr.a3('SWE');
    
    %get the final demand data for later weighting:
    QEXIO.Yglobal=sum(macro.yall.sec.all,2); % final demand
    QEXIO.Yswe=macro.yall.sec.all(:,strcmp('SWE',meta.countries)); % final demand
    
    
    cd(thisdir)
    %% Now we need to match the indicators to the stressors (and do any characterisation and unit conversion)
    clear tmp*
    for i=1:size(indic_hy.name,1)
        tmp_S=zeros(1,meta.Zdim);
        tmp_Fh=zeros(1,meta.Ydim);
        
        match_stress=strcmp(indic_hy.name(i),stressor_raw.Concordance.IndicatorNames);
        match_stress_orig_name=stressor_raw.Concordance.StressorNames(match_stress);
        match_stress_unit_conv=stressor_raw.Concordance.Characterisation(match_stress).*stressor_raw.Concordance.UnitConversion(match_stress);
        
        matched_stressors(i).name=indic_hy.name(i);
        matched_stressors(i).txt=[];
        matched_stressors(i).check_f=[];
        matched_stressors(i).unit_conv=[];
        disp(indic_hy.name(i))
        clear check_f
        for j=1:length(match_stress_orig_name)
            findx=find(strcmp(match_stress_orig_name(j),meta.labsF(:,1)));
            %             findx=find(strcmp(match_stress_orig_name(j),meta.labsF(:,2)) | strcmp(match_stress_orig_name(j),meta.labsF(:,1)));
            tmp_S=tmp_S+sum(stressors(findx,:),1)*match_stress_unit_conv(j);
            tmp_Fh=tmp_Fh+sum(stressors_hh(findx,:),1)*match_stress_unit_conv(j);
            matched_stressors(i).txt=[matched_stressors(i).txt;meta.labsF(findx,1)];
            if ~isempty(findx)
                matched_stressors(i).check_f=...
                    [matched_stressors(i).check_f;sum(stressors(findx,:),1)*match_stress_unit_conv(j)*IO.x];
                matched_stressors(i).unit_conv=...
                    [matched_stressors(i).unit_conv;match_stress_unit_conv(j)];
            end
        end
        
        
        %Intensities:
        QEXIO.S(:,i)=tmp_S;
        %Multipliers
        QEXIO.mult(:,i)=tmp_S*IO.L;
        %Multipliers by origin
        country_agg=kron(eye(49),ones(200,1));
%         QEXIO.mult_src(:,:,i)=tr.collapseZdim(bsxfun(@times,IO.L, QEXIO.S(:,i)),1); %same as:diag(QEXIO.S(:,i))*IO.L; % using mriotree
        QEXIO.mult_src(:,:,i)=country_agg'*(bsxfun(@times,IO.L, QEXIO.S(:,i))); %same as:diag(QEXIO.S(:,i))*IO.L; % using kronecker
        %HHld impact
        country_agg_y=kron(eye(49),ones(7,1));
%         tmp_Fh_cnt=tr.collapseYdim(tmp_Fh);
        tmp_Fh_cnt=(country_agg_y'*(tmp_Fh'))';
        QEXIO.F_hh(:,i)=tmp_Fh_cnt(sweden_index);
        %Sweden Production account:
        tmp_F=tmp_S.*IO.x';
%         tmp_F_cnt=tr.collapseZdim(tmp_F,2);
        tmp_F_cnt=(country_agg'*(tmp_F'))';
        QEXIO.prod(:,i)=tmp_F_cnt(sweden_index)+QEXIO.F_hh(:,i);
%         QEXIO.F(:,i)=tmp_F(tr.indexZ{sweden_index});
        QEXIO.F(:,i)=tmp_F(strcmp(meta.secLabsZ.Country,'SE'));
        %Sweden Consumption account:
%         tmp_Q_cnt=tr.collapseYdim(QEXIO.mult(:,i)'*IO.Y);
        tmp_Q_cnt=(country_agg_y'*(QEXIO.mult(:,i)'*IO.Y)')';
        QEXIO.cons(:,i)=tmp_Q_cnt(sweden_index);
        
        clear tmp_F tmp_F_cnt tmp_Q_cnt tmp_Fh_cnt tmp_S tmp_Fh
        
    end
    save(['QEXIO_',dataset_name,yrstr],'QEXIO','meta','sweden_index','matched_stressors')
end

%%
for i=1:size(indic_hy.name,1)
    if isempty(matched_stressors(i).txt)
        disp(['Error ' matched_stressors(i).name{1} , ' not matched'])
    end
end
