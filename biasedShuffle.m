function shuffledArray = biasedShuffle(originalArray, sigma1, sigma2, p)
    % This function does a biased shuffle of a given array
    % based on a given sigma hyperparameter
    N = length(originalArray);
    P = prctile(originalArray, p); % compute top p and bottom p percentiles


    shuffledArray = originalArray; % copy the original array

    % Normally distributed samples
    pd1 = makedist('Normal', 'mu', 0, 'sigma', sigma1);
    pd2 = makedist('Normal', 'mu', 0, 'sigma', sigma2);
    %samples = 
    for i = 1:N
        if or(shuffledArray(i) < P(1), shuffledArray(i) > P(2))
            samp = round(random(pd2, 1, 1)); % random number 
        else
            samp = round(random(pd1, 1, 1)); % get random number
        end
        j = i + samp;
        % check if j is one of the indices, if not set to i
        if or(j <= 0, j > N)
            j = i; % if out of range, j is i
        end
        % Swap value in index i and j
        temp = shuffledArray(i);
        shuffledArray(i) = shuffledArray(j);
        shuffledArray(j) = temp;
    end
end