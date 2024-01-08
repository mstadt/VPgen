function newArray = check_diffs(levORdsg, noOC, diff_lims)
    % This function ensures that the before OC and after OC differences
    % Inputs:
    %   levORdsg -- lev/dsg samples
    %   noOC -- noOC samples
    %   diff_lims -- min/max of the differences allowed

    % check the levORdsg - noOC values
    diff_OC = levORdsg - noOC;
    N = length(noOC);

    inds_outside_lim = find(diff_OC < diff_lims(1) | diff_OC > diff_lims(2));
    
    newArray = levORdsg;
    NUM_TRIALS = 0;
    while ~isempty(inds_outside_lim)
        NUM_TRIALS = NUM_TRIALS + 1;
        if NUM_TRIALS > 1e6
            fprintf('reached max trials! \n')
            fprintf('inds_outside_lim length: %i \n', length(inds_outside_lim))
            break
        end
       
    
        for i = 1:length(inds_outside_lim)
            ind = inds_outside_lim(i);

            noOC_val = noOC(ind);    
            OC_val = newArray(ind);

            OC_val_diffs = OC_val - noOC; % differences with OC val
            inds_within_range1 = find(OC_val_diffs >= diff_lims(1) & OC_val_diffs <= diff_lims(2));

            noOC_val_diffs = newArray - noOC_val; % diff with no OC val
            inds_within_range2 = find(noOC_val_diffs >= diff_lims(1) & noOC_val_diffs <= diff_lims(2));

            inds_within_range = intersect(inds_within_range1, inds_within_range2);

            if isempty(inds_within_range)
                %fprintf('NO POSSIBLE INDICES! \n')
                if length(inds_within_range1) < length(inds_within_range2)
                    inds_within_range = inds_within_range1;
                else
                    inds_within_range = inds_within_range2;
                end
            end


            rand_num = randi([1,length(inds_within_range)]);
            ind_j = inds_within_range(rand_num);
            
            % swap values
            temp = newArray(ind);
            newArray(ind) = newArray(ind_j);
            newArray(ind_j) = temp;
        end % for



        diff_OC = newArray - noOC;
        inds_outside_lim = find(diff_OC < diff_lims(1) | diff_OC > diff_lims(2));
        %fprintf('There are %i values outside limits \n', length(inds_outside_lim))
    end % while
    
    diff_OC = newArray - noOC;
    inds_outside_lim = find(diff_OC < diff_lims(1) | diff_OC > diff_lims(2));
    %if isempty(inds_outside_lim)
        %fprintf('all differences within limits! \n')
        %fprintf('number of trials: %i \n', NUM_TRIALS)
    %end
end