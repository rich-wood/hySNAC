        [a,chem_water_meta,c1]=xlsread('C:\Users\richardw\Box Sync\Projects\PRINCE_Swedish_EPA\WP4\InputData\20170530_PRTR_EMISSIONS_WATER_2013.xlsx');
        chem_water=a(4:end,4:7987+3);
        chem_water(isnan(chem_water))=0;
        for i=1:size(chem_water_meta,1)-3
            chem_water_names(i,1)={[chem_water_meta{i+3,2},'_water']};
            chem_water_names(i,3)=chem_water_meta(i+3,1);
            chem_water_names(i,2)=chem_water_meta(i+3,3);
        end        
        [a,chem_air_meta,c2]=xlsread('C:\Users\richardw\Box Sync\Projects\PRINCE_Swedish_EPA\WP4\InputData\20170530_PRTR_EMISSIONS_AIR_2013.xlsx');
        chem_air=a(4:end,4:7987+3);
        chem_air(isnan(chem_air))=0;
         for i=1:size(chem_air_meta,1)-3
            chem_air_names(i,1)={[chem_air_meta{i+3,2},'_air']};
            chem_air_names(i,3)=chem_air_meta(i+3,1);
            chem_air_names(i,2)=chem_air_meta(i+3,3);
        end 
        tox.val=[chem_water;chem_air];
        tox.names=[chem_water_names;chem_air_names];
        
        save('chem_tox_water_air','tox')
        
        char=importdata('C:\Users\richardw\Box Sync\Projects\PRINCE_Swedish_EPA\WP4\characterisation\Compiled_UseTox2.xlsx');
        %%
        
        air_result=chem_air'*char.data.AIR(1:45,:);
        water_result=chem_water'*char.data.Water(:,:);
        
        tox.val=(air_result+water_result)';
        tox.names=[char.textdata.Water(1,3:4)',char.textdata.Water(1,3:4)',char.textdata.Water(1,3:4)']
        
        save('chem_tox_water_air_char','tox')
      
        
        swe_indx=(strcmp('SE',chem_air_meta(2,4:end)));
        
        air_result_swe=sum(air_result(swe_indx,:))
        water_result_swe=sum(water_result(swe_indx,:))
        air_result_global=sum(air_result(:,:))
        water_result_global=sum(water_result(:,:))
        
        swe_result=air_result_swe+water_result_swe
        
        
        