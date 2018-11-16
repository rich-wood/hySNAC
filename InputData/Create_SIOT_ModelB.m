function [T,A,Tdom,Timp,x,Y, Adom,Aimp] = Create_SIOT_ModelB(sup,usedom,useimp)
% #constructs
% this script allocates noncompetitive imports as an additional input row
% (not in the intermediate industry classification)
% this becuase it is impossible to substitute a noncompetitive import
% domestically.


indout = sum(sup,1);
n_p=size(sup,1);
n_ind=size(sup,2);
indout(indout<1e-5)=1;
prodout_dom = sum(sup,2);
prodout_dom(prodout_dom<1e-5)=1;

if isempty(useimp)
    useimp=usedom*0;
end

usetot=usedom;
usetot(1:n_p,:)=usetot(1:n_p,:)+useimp(1:n_p,:);

make=transpose(sup);


%(Eurostat Model B)
% industry technology assumption
% Product-by-Product input-output table - no negatives guaranteed
%     disp('negatives in A, using alternate technology assumption')
%     D=inv(diag(indout))*(make); same as next line in matlab notation:
D=diag(indout)\(make);
Tdom=usedom(:,1:n_ind)*D;
Timp=useimp(1:n_p,1:n_ind)*D;
Timp(sum(sup,2)==0,:)=0;

T=Tdom;
T(1:n_p,:)=T(1:n_p,1:n_p)+Timp(:,1:n_p);

A=T(:,1:n_p)/(diag(prodout_dom));
Adom=Tdom(:,1:n_p)/(diag(prodout_dom));
Aimp=Timp(:,1:n_p)/(diag(prodout_dom));


Y=usedom(1:n_p,n_ind+1:end)+useimp(1:n_p,n_ind+1:end);

prodout_dom = sum(sup,2);
x=prodout_dom;
