%% 
clc; close all;
load Subject1.mat
load Subject2.mat

Sub1 =transpose(table2array(subject1));
Sub2 =transpose(table2array(subject2));

%Sub1_1 = Sub1(1:19,:);
%sub2_1 = Sub2(1:19,:);

save Sub1.mat;
save Sub2.mat;


% function epochMatrix = epoch(eegData)
% 
% data = EEG.data;
% for j = 1 : 19
%     for i = 1 : 120
%         epochMatrix(j,:,i) = data(j,(i*2000-599):(i*2000));    
%     end        
% end
% 
% 
% end‏
% ‬
% 
% UnsampledEpoch = EEG.data;
% subSampledEpoch2 = subSample(UnsampledEpoch);
% outputEpoch2.epoch = subSampledEpoch2;
% outputEpoch2.oder = normal(1).odor(1:118,1);
% outputEpoch2.noisy(1,1) = 1;
% outputEpoch2.noisy(1,2) = 2;
% outputEpoch2.noisy(1,3) = 3;
% save outputEpoch2.mat
% function subSampledEpoch = subSample(UnsampledEpoch)
% This function prepare subsamples of the epoch

% subSampledEpoch2(1,:,:) = UnsampledEpoch(1,:,:);
% subSampledEpoch2(2,:,:) = UnsampledEpoch(5,:,:);
% subSampledEpoch2(3,:,:) = UnsampledEpoch(10,:,:);
% subSampledEpoch2(4,:,:) = UnsampledEpoch(15,:,:);‏‏