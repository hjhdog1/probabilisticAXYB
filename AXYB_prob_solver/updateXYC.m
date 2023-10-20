function [X,Y,C] = updateXYC(X, Y, C, GwX, GwY, GwC, GqX, GqY, GqC, step_R, step_p, nData)

epsilon = 1e-6;

% X = X * [LargeSO3(step_R/nData*GwX), step_p/nData*GqX'; 0,0,0,1];
% Y = Y * [LargeSO3(step_R/nData*GwY), step_p/nData*GqY'; 0,0,0,1];
X = X * [LargeSO3(step_R/nData*GwX/(norm(GwX)+epsilon)), step_p/nData*GqX'/(norm(GqX)+epsilon); 0,0,0,1];
Y = Y * [LargeSO3(step_R/nData*GwY/(norm(GwY)+epsilon)), step_p/nData*GqY'/(norm(GqY)+epsilon); 0,0,0,1];

n = size(C,3);
for i = 1:n
% %     C(:,:,i) = C(:,:,i) * [LargeSO3(step_R*GwC(i,:)), step_p*GqC(i,:)'; 0,0,0,1];
%     C{i} = C{i} * [LargeSO3(step_R*GwC(i,:)), step_p*GqC(i,:)'; 0,0,0,1];
    C{i} = C{i} * [LargeSO3(step_R*GwC(i,:)/(norm(GwC(i,:))+epsilon)), step_p*GqC(i,:)'/(norm(GqC(i,:))+epsilon)'; 0,0,0,1];
end


% 
% X = X * [LargeSO3(step_R/nData*GwX), step_p/nData*GqX'; 0,0,0,1];
% Y = Y * [LargeSO3(step_R/nData*GwY), step_p/nData*GqY'; 0,0,0,1];
% 
% n = size(C,3);
% for i = 1:n
%     C{i} = C{i} * [LargeSO3(step_R*GwC(i,:)), step_p*GqC(i,:)'; 0,0,0,1];
% end
% 
end

