function [acc,feat_weights] = do_my_classification(vectors_train,labels_train,vectors_test,labels_test,STUDY)
%
% Performs multivariate pattern classification/regression using input
% vectors of training/test data and condition labels. 
%
% This script interacts with LIBSVM toolbox (Chang & Lin) to do the classfication / regression
% see: https://www.csie.ntu.edu.tw/~cjlin/libsvm/
% Chang CC, Lin CJ (2011). LIBSVM : a library for support vector machines. ACM TIST, 2(3):27,
%
% 
% Inputs:
%
%   vectors_train   data vectors that make up the training dataset
%   labels_train    condition labels for the training dataset
%   vectors_test    data vectors that make up the test dataset
%   labels_test     condition labels for the test dataset
%   STUDY           structure containing participant dataset information and 
%                   multivariate classification/regression settings.
% 
%
% Outputs:
%
%   acc             classifier accuracy for classification of the test dataset.
%   feat_weights    feature weights from the multivariate pattern
%                   classification.
%
%
% Copyright (c) 2013-2016 Stefan Bode and contributors
% 
% This file is part of DDTBOX.
%
% DDTBOX is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% define samples and labels for training

Samples = vectors_train;
Labels = labels_train;


%% training 
%__________________________________________________________________________
if sum(STUDY.analysis_mode == [1 3 4]) %  libsvm
    model = svmtrain(Labels,Samples,STUDY.backend_flags.all_flags);
elseif sum(STUDY.analysis_mode == [2]) %  liblinear
   model = train(Labels,sparse(Samples),STUDY.backend_flags.all_flags);
end
%__________________________________________________________________________    


%% define samples and labels for testing
Samples = vectors_test;
Labels = labels_test;


%% prediction
%__________________________________________________________________________

if sum(STUDY.analysis_mode == [1 3]) % libsvm
    [predicted_label, accuracy, decision_values] = svmpredict(Labels, Samples, model); 
elseif sum(STUDY.analysis_mode == [2]) % liblinear
    [predicted_label, accuracy, decision_values] = predict(Labels, sparse(Samples), model); 
end

if STUDY.analysis_mode == 1 % SVM classification with libsvm
    % calculating feature weights
    w = model.SVs' * model.sv_coef;
    b = -model.rho;

    feat_weights = zeros(size(w,1),3);
    
    if STUDY.feat_weights_mode == 1 
        feat_weights(:,1) = 1:(size(w,1));
        feat_weights(:,2) = w;
        feat_weights(:,3) = abs(w);
    end

    % extracting accuracy for 2-classes 
    if STUDY.nconds == 2

        acc = accuracy(1);

    % extracting accuracy for N-classes
    elseif STUDY.nconds > 2

        classes = 1:STUDY.nconds;
        pairs = nchoosek(classes,2);

        for cl = classes
            wt(:,cl) = (pairs(:,1) == cl) - (pairs(:,2) == cl);
        end

        votes = decision_values * wt;
        [maxvote,winvote] = max(votes');
        classcorrectnes = (classes == winvote) * 100;

        acc = mean(classcorrectnes);

    end % if nclass
    
elseif STUDY.analysis_mode == 2 % SVM classification with liblinear
    
    % calculating feature weights
    %w = model.SVs' * model.sv_coef;
    %b = -model.rho;
    w = model.w';

    feat_weights = zeros(size(w,1),3);
    
    if STUDY.feat_weights_mode == 1 
        feat_weights(:,1) = 1:(size(w,1));
        feat_weights(:,2) = w;
        feat_weights(:,3) = abs(w);
    end

    % extracting accuracy for 2-classes 
    if STUDY.nconds == 2

        acc=accuracy(1);
    
    % extracting accuracy for N-classes
    elseif STUDY.nconds > 2

        classes = 1:STUDY.nconds;
        pairs = nchoosek(classes,2);

        for cl = classes
            wt(:,cl) = (pairs(:,1)==cl) - (pairs(:,2)==cl);
        end

        votes = decision_values*wt;
        [maxvote,winvote] = max(votes');
        classcorrectnes = (classes==winvote)*100;

        acc = mean(classcorrectnes);

    end % if nclass
       
elseif STUDY.analysis_mode == 3 % SVR
    
    % correlating the predicted label with the test labels
    c_sample = corrcoef(predicted_label,Labels);
    avecorrectness = c_sample(1,2);
    
    % convert into Fisher-Z
    correctness_z = 1/2 * log((1 + avecorrectness) ./ (1 - avecorrectness));
    avecorrectness = mean(correctness_z);
    
    % optional: transform back
    % avecorrectness = (exp(2 * avecorrectness) - 1) ./ (exp(2 * avecorrectness) + 1);
    
    acc = avecorrectness;
    
    % calculating feature weights
    w = model.SVs' * model.sv_coef;
    feat_weights = zeros(size(w,1),3);  
    
    if STUDY.feat_weights_mode == 1
        feat_weights(:,1) = 1:(size(w,1));
        feat_weights(:,2) = w;
        feat_weights(:,3) = abs(w);
    end
    
end % if analysis_mode
%__________________________________________________________________________
