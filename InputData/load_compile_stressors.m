
function [F,Fhh]=load_compile_stressors(indic,yr,QEXIO,SWE,C_prod)


yrstr=num2str(yr);



for i=1:length(indic.name)
    if strcmp(indic.src(i),'EXIOBASE')
        Fbig=QEXIO.F(:,i); %Environmental extensions to industry usage
        Fhh(i,:)=QEXIO.F_hh(i)'; %Environmental extensions to final usage
        F(i,:)=Fbig'*C_prod'; %Aggregate to swedish classification
    elseif strcmp(indic.src(i),'SWE')
        swe_indx=nnz(strcmp(indic.src(1:i),"SWE"));
        F(i,:)=SWE.F(swe_indx,:); %Environmental extensions to industry usage
        Fhh(i,:)=SWE.Fhh(swe_indx,:)'; %Environmental extensions to final usage
    end
end

end

