function [train,test,ntrain,ntest] = wdbc(datafile,dataDim,fracTest,reord)
% syntax: [train,test,ntrain,ntest] = wdbcData(datafile,dataDim,fracTest,reord)
% extract data from the database
% here, "datafile" should be a string, eg 'wdbc.data'
%       "dataDim" is a scalar, eg 30
%       "fracTest" is a scalar strictly between 0 and 1 indicating
%       what fraction of data should go in the test set (the
%       remaining data goes in the training set)
%       "reord" indicates whether the data should be reordered
%       before selecting the test/training sets. Value of "0"
%       indicates no reordering, "1" indicates random reordering
%
% on return, ntrain and ntest indicate the number of rows in the
% training and testing sets, respectively. train and test contain the
% data arrays. Elements in first column of these matrices are either 0
% (indicating a benign sample) or 1 (indicating a malignant
% sample). Elements in the remaining columns, starting with column
% 2, contain the features.


% first check input data
if (nargin < 3)
    error('three or four input arguments are required for wdbcData');
end
if (~ischar(datafile))
    error('first argument must be a string');
end
if (isnumeric(dataDim)&max(size(dataDim))==1)
    if (isempty(dataDim))
        error('second argument must be a scalar');
    elseif (dataDim <= 0)
        error('second argument must be a positive integer');
    end
else
    error('second argument must be a positive number');
end
if(isnumeric(fracTest))
%     if(fracTest<=0 | fracTest>=1)
%         error('third argument must be a number strictly between 0 and 1');
%     end
else
    error('third argument must be a number');
end

if nargin==4
    if(~isnumeric(reord))
        error('fourth argument should be numeric, either 0 or 1');
    end
else
    % default is no reordering
    reord = 0;
end
reord = 0;



fp = fopen(datafile,'r');
samples = 0;
train = zeros(0,dataDim+1);
[id,count] = fscanf(fp,'%d,',1);
while (count == 1)
    samples = samples + 1;
    type = fscanf(fp,'%1s',1);
    if type=='B'
        type =0;
    elseif type=='M'
        type =1;
    else
        type=-1;
        fprintf(' invalid type found in data file %s\n', datafile);
    end
    vec = fscanf(fp,',%e',dataDim);
    train = [train; type vec'];
    [id,count] = fscanf(fp,'%d,',1);
end
if (samples < 569)
    error('Not enough samples');
end

% reorder the rows of "train", if requested
if ~(reord==0)
    p = randperm(samples);
    train(p,:) = train;
end

ntest = round(fracTest * samples);
if ntest < 1
    ntest=1;
elseif ntest >= samples
    ntest = samples-1;
end
ntrain = samples - ntest;

% test = train(1:ntest,:);
% train = train(ntest+1:samples,:);
% modified 4/7/03
test  = train(ntrain+1:samples,:);
train = train(1:ntrain,:);
end
