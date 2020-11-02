function [jfound,idx_1,idx_2] = fn_jumpfind(input_vect, jumpsize)
    %finds first idx_1 and idx_2 where
    %input_vect(idx_2)-input_vect(idx_1)>= jumpsize
    jfound = 0;
    idx_1 = 1;
    idx_2 = length(input_vect);
    final_jfound = 0;
    
    idx_pivot = ceil(idx_2/2);
    idx_t1 = 1;
    idx_t2 = 1;
    
    if(input_vect(idx_2)-input_vect(idx_1) >= jumpsize)
        jfound = 1;
    end
    
    while(final_jfound ~= 1)
        if(jfound == 1)
            % search in this one
            if(idx_2-idx_1 == 1)
                % big Yes 
                %that's it, we found it
                final_jfound = 1;
            else
                % well yes but not really
                %find new idxes, gotta go deeper
                %search left half
                [jfound,idx_t1,idx_t2] = fn_jumpfind(input_vect(idx_1:idx_pivot), jumpsize);
                if(jfound ~= 1)
                    %no jfound, look at right of pivot instead
                    [jfound,idx_t1,idx_t2] = fn_jumpfind(input_vect(idx_pivot:idx_2), jumpsize);
                    if(jfound ~= 1)
                        %no jfound again so this vect has no jump at all
                        final_jfound = 1;
                    else
                        %jfound in right side, update idxes
                        idx_1 = idx_t1 + idx_pivot - 1;
                        idx_2 = idx_t2 + idx_pivot - 1;
                        if(idx_2-idx_1 == 1)
                            %that's it, we found it
                            final_jfound = 1;
                        end
                    end
                else
                    %jfound on left side, update idxes
                    idx_1 = idx_t1 + idx_1 - 1;
                    idx_2 = idx_t2 + idx_1 - 1;
                    if(idx_2-idx_1 == 1)
                        %that's it, we found it
                        final_jfound = 1;
                    end
                end
            end
        else
           %no jfound in this one, move up and update indexes
           final_jfound = 1;
        end
    end
end



