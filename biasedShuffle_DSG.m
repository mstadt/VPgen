function shuffledArray = biasedShuffle_DSG(originalArray, sigma1, sigma2, p, diff_lev, noOC)
    % Additional assumption: DSG is biased based on Lev effects 
    % Inputs: diff_lev -- lev differences, noOC -- noOC array

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
        % Check if direction matches lev direction
        % If not, find a pair that matches lev direction




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
        % Check if values are in same direction as lev change
        current_val = shuffledArray(i);
        new_val = shuffledArray(j);

        % current difference with current value
        diff_i = sign(current_val - noOC(i));
        diff_i_lev = sign(diff_lev(i));
        if diff_i ~= diff_i_lev
            flag = 1;
        end
        
        % difference at current index with new value
        new_diff_i = sign(new_val - noOC(i));

        % difference at new index
        new_diff_j = sign(current_val - noOC(j));

        if and(new_diff_i == sign(diff_lev(i)), new_diff_j == sign(diff_lev(j)))
            % swap value in index i and j if the new directions
            %  match the lev directions
            temp = shuffledArray(i);
            shuffledArray(i) = shuffledArray(j);
            shuffledArray(j) = temp;
        end
        

    end
end