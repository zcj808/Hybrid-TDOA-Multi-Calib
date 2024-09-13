function [mat] = skew(vec)
    mat=zeros(3);
    mat(1,2)=-vec(3);
    mat(1,3)=vec(2);
    mat(2,3)=-vec(1);
    mat=mat-mat';
end

