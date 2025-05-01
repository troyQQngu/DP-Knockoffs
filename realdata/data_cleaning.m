%% 1. Load raw data
data = readtable('diabetic_data.csv');

% 1) Ensure age is a string array
if iscell(data.age)
    ageStr = string(data.age);
elseif iscategorical(data.age)
    ageStr = string(categories(data.age));
    % Map back into a full vector
    ageStr = ageStr(data.age);  
else
    ageStr = data.age;  % already string array
end

% 2) Precompute mid‐points for each unique category
uniqCats = unique(ageStr);
midPoints = zeros(size(uniqCats));

for i = 1:numel(uniqCats)
    s = uniqCats(i);
    % try to match “[low,high)”
    tokens = regexp(s, '\[(\d+)-\s*(\d+)\)', 'tokens');
    if ~isempty(tokens)
        low  = str2double(tokens{1}{1});
        high = str2double(tokens{1}{2});
        mid  = (low + high) / 2;
    else
        % fallback for “>90” or other formats
        num = sscanf(s, '>%d');
        if ~isempty(num)
            mid = num + 5;    % e.g. “>90” → 95
        else
            mid = NaN;        % or some default
        end
    end
    midPoints(i) = mid;
end

% 3) Build a lookup map from category→midpoint
ageMap = containers.Map( uniqCats, midPoints );

% 4) Apply it to the full age vector
n = height(data);
ageNum = zeros(n,1);
for j = 1:n
    ageNum(j) = ageMap(ageStr(j));
end

% 5) Replace in your table
data.age_numeric = ageNum;
data.age          = [];    % drop the old categorical/string

% 6) (Optional) Standardize the numeric “age_numeric”
data.age_numeric = (data.age_numeric - mean(data.age_numeric, 'omitnan')) ...
                                   / std(data.age_numeric,  'omitnan');


%% 2. Drop obvious non-predictors
data(:,{'patient_nbr','encounter_id','admission_type_id','discharge_disposition_id','admission_source_id'}) = [];    % IDs
data.weight = [];                                % > 97% missing

%% 3. Drop any column with >90% missing values
missingPct = sum(ismissing(data)) / height(data);
colsHighMissing = data.Properties.VariableNames(missingPct > 0.90);
data(:, colsHighMissing) = [];

%% 4. Separate response
y = data.time_in_hospital;       % 1–14 days
data.time_in_hospital = [];

%% 5. Identify categorical vs. numeric
varTypes = varfun(@class, data, 'OutputFormat', 'cell');
isCat    = strcmp(varTypes,'cell') | strcmp(varTypes,'char') | strcmp(varTypes,'string');
catVars  = data.Properties.VariableNames(isCat);
numVars  = setdiff(data.Properties.VariableNames, catVars);

%% 6. Fill missing in categoricals and convert to categorical
for v = catVars
    col = data.(v{1});
    col = fillmissing(col, 'constant', 'Unknown');
    data.(v{1}) = categorical(col);
end

%% 7. Fill & standardize numeric predictors
for k = 1:numel(numVars)
    name = numVars{k};
    col  = data.(name);

    % 1) Compute median ignoring NaNs
    m = median(col, 'omitnan');

    % 2) Replace any NaNs with the median
    col = fillmissing(col, 'constant', m);

    % 3) Now z-score (omitnan no longer needed since no more NaNs)
    data.(name) = (col - mean(col)) / std(col);
end

%% 8. Expand categoricals into dummy variables *and* standardize them
for v = catVars
    % original categorical column
    C = data.(v{1});
    
    % 1) one–hot encode
    D = dummyvar(C);              % size [n×#levels]
    
    % 2) drop the first level (reference)
    D(:,1) = [];                  % now [n×(#levels−1)]
    
    % 3) standardize each dummy column
    %    (zscore by column)
    D = zscore(D);                % subtract mean, divide by std
    
    % 4) assign into the table
    lvl = categories(C);
    newNames = strcat(v{1}, '_', lvl(2:end));
    for j = 1:numel(newNames)
        data.(newNames{j}) = D(:,j);
    end
    
    % 5) remove the original categorical
    data.(v{1}) = [];
end


%% 9. Remove any zero-variance predictors
stdVals = varfun(@std, data, 'OutputFormat', 'uniform');
zeroVar  = data.Properties.VariableNames(stdVals == 0);
data(:, zeroVar) = [];

%% 10. Build final matrices and save
X = table2array(data);
save('diabetic_clean.mat', 'X', 'y','data');

