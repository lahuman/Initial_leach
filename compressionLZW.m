function [packLen] = compressionLZW (dataIn)

originalText = [dataIn]; % Makes the "originalText" a char class variable from cell class.

%originalText = cell2mat(originalText);
initialDictionary = unique(originalText); % Chooses the unique characters of the input text.
fprintf('Original Text:%s\n', originalText); % Displays the original text.
fprintf('Initial Dictionary:%s\n', initialDictionary); % Displays the initial dictionary.
for c = 1 : length(initialDictionary) % This loop counts the number of iterations of each unique character of the input text.
	countForThisChar = sum(originalText == initialDictionary(c));
	fprintf(2,'Number of iterations of character %s : %d\n', initialDictionary(c), countForThisChar); % Displays number of iteratios in red color.
end


Dictionary = cell(length(initialDictionary),1);
for j=1:length(initialDictionary)
    Dictionary(j) = {initialDictionary(j)};
end

P = '';
P_code = -1;
k = 0;
idicmaxLen = length(initialDictionary);

for i=1:length(originalText)
    Q = strcat(P,originalText(i));
    for j=1:idicmaxLen
        Q_code = 0;
        if (strcmp(Q, Dictionary(j)) == 1)
            Q_code = j;
            break;
        end
    end
    if (Q_code > 0)
        P = Q;
        P_code = Q_code;
    else
        k = k+1;
        output(k) = P_code;
        idicmaxLen = idicmaxLen + 1;
        Dictionary(idicmaxLen) = {Q};
        P = originalText(i);
        P_code = strfind(initialDictionary, P);
    end
end
k = k+1;
output(k) = P_code;
display(output);

[packLen, n]=  size(dec2bin((sprintf('%s', output).'))-'0');
% str_x = num2str(output);
% str_x(isspace(str_x)) = ''

% packLen = dec2bin(output);

end