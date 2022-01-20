function [ARR_ROUTH,result]=JustifyRough(poly)

%part1 calculate the routh arARR_ROUTHy
%
%   The first part will calculate the symbolic Routh array for 
%   characteristic polynomial. The following special cases are considered:

%   1) zero first elements and 2) rows of zeros. All zero first 
%   elements are replaced with the variable EPSILON
%   which can be later substituted with positive small
%   numbers[in this code,small number is 1/100 * min(poli)]. When a  
%   row of zeros is found, the auxiliary polynomial is used.
%
%   Examples:
%
%   1) Routh array for s^3+2*s^2+3*s+1
%
%       >>ARR_ROUTH=routh([1 2 3 1])
%
%   2) Routh arARR_ROUTHy for s^3+a*s^2+b*s+c
%   
%       >>syms a b c;
%       >>ARR_ROUTH=routh([1 a b c]);
%
% 
%
flag = 0;    %to help judging the result,if row of zeros happen,set flag equals to 1
epsilon = min(poly)/100; %select one percent of the minimum of poli as epsilon
% if(nargin<2),
%   fprintf('\nError: Not enough input arguments given.');
%   return
% end
if all(poly<0)
    poly = -1 * poly;       %tARR_ROUTHnsform poli to all positive
end

dim=size(poly);     %get size of poly       

coeff=dim(2);               %get number of coefficients
ARR_ROUTH=sym(zeros(coeff,ceil(coeff/2))); %initialize symbolic Routh arARR_ROUTHy 

for i=1:coeff,
    ARR_ROUTH(2-rem(i,2),ceil(i/2))=poly(i); %assemble 1st and 2nd rows
end

rows=coeff-2;       %number of rows that need determinants
index=zeros(rows,1);    %initialize columns-per-row index vector

for i=1:rows,
    index(rows-i+1)=ceil(i/2); %form index vector from bottom to top
end

for i=3:coeff,              %go from 3rd row to last
    if(all(ARR_ROUTH(i-1,:)==0)),      %row of zeros
            flag = 1;           %to help judging result
            fprintf('\nSpecial Case: Zero value in row detected.\n');
            a=coeff-i+2;        %order of auxiliary equation
            b=ceil(a/2)-rem(a,2)+1; %number of auxiliary coefficients
            temp1=ARR_ROUTH(i-2,1:b);  %get auxiliary polynomial
            temp2=a:-2:0;       %auxiliry polynomial powers
            ARR_ROUTH(i-1,1:b)=temp1.*temp2;   %derivative of auxiliary
    elseif(ARR_ROUTH(i-1,1)==0),       %first element in row is zero
            fprintf('\nSpecial Case: The first element is zero.\n');
            ARR_ROUTH(i-1,1)=epsilon;  %replace by epsilon
    end
                %compute the Routh arARR_ROUTHy elements
    for j=1:index(i-2), 
        ARR_ROUTH(i,j)=-det([ARR_ROUTH(i-2,1) ARR_ROUTH(i-2,j+1);ARR_ROUTH(i-1,1) ARR_ROUTH(i-1,j+1)])/ARR_ROUTH(i-1,1);
    end
end

%The second part: Based on the Routh array, justify the stablity of the
%whole system
first_colomn_ARR_ROUTH = ARR_ROUTH(:,1);
if min(first_colomn_ARR_ROUTH)<0
    result = 'The system is unstable';
elseif all(first_colomn_ARR_ROUTH>0)      %first_colomn_ARR_ROUTH has no zero
    if flag == 1
        result = 'The stablity can not be judged simply by Routh Criterion';
    elseif flag == 0
        result = 'The system is stable';
    end
end