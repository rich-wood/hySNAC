function prep_swe_data(yr)

yrstr=num2str(yr);

%%     DATA READ AND SET UP

%% load all the data, set up meta, and create indices for Z and y matrices for SWE_srcden only
SWE_src=importdata(['InputData\IO_SE_08_14_170607.xlsx']); % just a variable of all data and meta.
SWE.Ud=SWE_src.data.(['x',yrstr])(3+1:3+59,2+1:2+59); % the domestic use table
SWE.Um=SWE_src.data.(['x',yrstr])(85+1:85+59,2+1:2+59); % the import use table
SWE.V=SWE_src.data.(['x',yrstr])(168+1:168+59,2+1:2+59); % the supply table
SWE.Yd=SWE_src.data.(['x',yrstr])(3+1:3+59,68:75); % the domestic final demand
SWE.Ym=SWE_src.data.(['x',yrstr])(85+1:85+59,68:75); % the import final demand
SWE.PI=SWE_src.data.(['x',yrstr])(68+1:68+59,2+1:2+59); % primary inputs
SWE.trade_shares=SWE_src.data.(['x',yrstr])(239:297,2+1:204); % the import use table
SWE.M=sum(SWE.Um,2)+sum(SWE.Ym,2);
%NOTE - environmental data read in lower down

SWE.cnt_names=SWE_src.textdata.x2008([235,237],3:end);
SWE.prod_names=SWE_src.textdata.x2008(168:227,2);
SWE.Fd_names=SWE_src.textdata.x2008(2,68:75);

% Breakdown of Imports by country of origin:
SWE.trade_by_country_of_origin_by_product=SWE.trade_shares'*diag(SWE.M);

% Create symmetric IO tables (uses existing scripts):
[SWE.T,SWE.A,SWE.Tdom,SWE.Timp,SWE.x,SWE.Y,SWE.Adom,SWE.Aimp] = Create_SIOT_ModelB(SWE.V,SWE.Ud,SWE.Um);
% Create coefficient versions


% Now create Ldom (inverse of diagonalised Adom)
% for a single region model, this is just, e.g. 59x59, and not a full MRIO
% layout
SWE.Ldom=inv(eye(size(SWE.Adom))-SWE.Adom);
%read SWE_srcdish F data from excel
SWE.F=SWE_src.data.(['x',yrstr])([150:163,71],2+1:2+59); %Environmental extensions to industry usage
SWE.Fhh=SWE_src.data.(['x',yrstr])([150:163,71],64); %Environmental extensions to final usage (&value added )
SWE.Fhh(isnan(SWE.Fhh))=0;
SWE.Fnames=SWE_src.textdata.(['x',yrstr])([149:162,70],1:2);

save(['SWE',yrstr],'SWE')

end
