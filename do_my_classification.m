function [acc,feat_weights, feat_weights_corrected] = do_my_classification(vectors_train,labels_train,vectors_test,labels_test,STUDY)
%__________________________________________________________________________
% DDTBOX script written by Stefan Bode 01/03/2013
%
% The toolbox was written with contributions from:
% Daniel Bennett, Daniel Feuerriegel, Phillip Alday
%
% The author further acknowledges helpful conceptual input/work from: 
% Jutta Stahl, Simon Lilburn, Philip L. Smith, Elaine Corbett, Carsten Murawski, 
% Carsten Bogler, John-Dylan Haynes
%__________________________________________________________________________
%
% This script interacts with LIBSVM toolbox (Chang & Lin) to do the classfication / regression
% see: https://www.csie.ntu.edu.tw/~cjlin/libsvm/
% Chang CC, Lin CJ (2011). LIBSVM : a library for support vector machines. ACM TIST, 2(3):27,
%
%
%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

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
    feat_weights_corrected = zeros(size(w,1),3);
    
    if STUDY.feat_weights_mode == 1 
        % uncorrected feature weights
        feat_weights(:,1) = 1:(size(w,1));
        feat_weights(:,2) = w;
        feat_weights(:,3) = abs(w);
        
        % corrected feature weights according to Haufe et al. (2014) method
        feat_weights_corrected(:,1) = 1:(size(w,1));
        vectors_train_temp = (vectors_train - repmat(mean(vectors_train), size(vectors_train, 1), 1)) ./ sqrt(size(vectors_train, 1) - 1);
        feat_weights_corrected(:,2) = vectors_train_temp' * (vectors_train_temp * feat_weights(:,2));
        clear vectors_train_temp;
        feat_weights_corrected(:,3) = abs(feat_weights_corrected(:,2));
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
    feat_weights_corrected = zeros(size(w,1),3);
    
    if STUDY.feat_weights_mode == 1 
        % uncorrected feature weights
        feat_weights(:,1) = 1:(size(w,1));
        feat_weights(:,2) = w;
        feat_weights(:,3) = abs(w);
        
        % corrected feature weights according to Haufe et al. (2014) method
        feat_weights_corrected(:,1) = 1:(size(w,1));
        vectors_train_temp = (vectors_train - repmat(mean(vectors_train), size(vectors_train, 1), 1)) ./ sqrt(size(vectors_train, 1) - 1);
        feat_weights_corrected(:,2) = vectors_train_temp' * (vectors_train_temp * feat_weights(:,2));
        clear vectors_train_temp;
        feat_weights_corrected(:,3) = abs(feat_weights_corrected(:,2));
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
    feat_weights_corrected = zeros(size(w,1),3);
    
    if STUDY.feat_weights_mode == 1
        % uncorrected feature weights
        feat_weights(:,1) = 1:(size(w,1));
        feat_weights(:,2) = w;
        feat_weights(:,3) = abs(w);
        
        % corrected feature weights according to Haufe et al. (2014) method
        feat_weights_corrected(:,1) = 1:(size(w,1));
        vectors_train_temp = (vectors_train - repmat(mean(vectors_train), size(vectors_train, 1), 1)) ./ sqrt(size(vectors_train, 1) - 1);
        feat_weights_corrected(:,2) = vectors_train_temp' * (vectors_train_temp * feat_weights(:,2));
        clear vectors_train_temp;
        feat_weights_corrected(:,3) = abs(feat_weights_corrected(:,2));
    end
    
end % if analysis_mode
%__________________________________________________________________________
