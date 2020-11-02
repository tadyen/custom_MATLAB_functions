function [out1, out2] = fn_linterpolate(in1, in2, in3)
    % by Alex 'tadyen' Wong, 24-10-2019
    % version 2.2 : 4-12-2019
    % Modes:
    % 1) [Out_X_new, ~] = fn( X(t), itpl_upsample_factor ); 
    %       args: 2
    % 2) [Out_X_new, Out_Y_new] = fn( X(t), Y(t), intp_upsample_factor );
    %       args: 3
    % 3) [Out_Y_new, ~] = fn( X_old(t), X_new(t), Y(t) );
    %       args: 3 , Native mode in internal function
    %
    % Resamples 1D-arrays by Linearly interpolating them. Arr_Y is linearly
    % interpolated to closeness between respective arr_X points.
    % Inputting only 1 array results in a linear interpolation against
    % itself, running the assumption that each point is evenly spaced in
    % it's parametric axis.
    

    mode = 0;
    
    switch nargin
        case 3
            
            % Selecting the mode
            if(max(size(in3)~=[1,1])==1)
                %Mode 3: [Out_Y_new, ~] = fn( X_old(t), X_new(t), Y(t) );
                if(min(size(in1))~= 1 || min(size(in2))~=1 || min(size(in3))~=1)
                    error('Input has incorrect size');
                else
                    if(length(in1)~=length(in3))
                        error('Input has incorrect size');
                    end
                end
                mode = 3;
            else
                %Mode 2) [Out_X_new, Out_Y_new] = fn( X(t), Y(t), intp_upsample_factor );
                if(min(size(in3))~= 1)
                    error('Input has incorrect size');
                else
                    if(min(size(in1))~= 1 || min(size(in2))~=1)
                        error('Input has incorrect size');
                    else
                        if(length(in1)~=length(in2))
                            error('Input has incorrect size');
                        end
                    end
                    mode = 2;
                end
            end
            
            
            % Performing action based on mode
            switch mode
                case 2
                    %2) [Out_X_new, Out_Y_new] = fn( X(t), Y(t), intp_upsample_factor );
                    out1 = fn_linterpolate(in1, in3);
                    out2 = fn_linterpolate_internal(in1, out1, in2);
                    
                case 3
                    %Mode 3: [Out_Y_new, ~] = fn( X_old(t), X_new(t), Y(t) );
                    out1 = fn_linterpolate_internal(in1, in2, in3);
                    out2 = [];
                otherwise
                    error('Houston weve got a Jim its dead');
            end
            
            
            
        case 2
            % 1) Out_X_new = fn( X(t), intp_upsample_factor ); 
            %    args: 2
            
            if(max(size(in2)~=[1,1])==1)
                error('Incorrect input format');
            else
                if(in2 <= 0)
                    error('Upsample factor must be greater than 0');
                end
                itpl_upsample_factor = in2;
            end
            in_Y_old = in1;
            
            idxlen_1 = length(in1);
            in_X_old = 1:idxlen_1;
            
            idxlen_2 = ceil(idxlen_1*itpl_upsample_factor)-1;
            in_X_new = linspace(1,idxlen_1,idxlen_2);
            
            out1 = fn_linterpolate_internal(in_X_old, in_X_new, in_Y_old);
        case 1
            error('not enough input args');
        otherwise
            error('Houston weve got a Jim its dead');
    end
    
    
    
    
    
    
    % ==================================================================
    % ===                 Internal Function                          ===
    % ==================================================================
    function Y_new = fn_linterpolate_internal(X_old, X_new, Y_old)
        %  Y_new = fn( X_old(t), X_new(t), Y(t) );
        intl_idxlen_1 = length(X_old);
        intl_idxlen_2 = length(X_new);
        Y_new = zeros(intl_idxlen_2, 1);
        
        if(min(X_old) <= min(X_new))
            in_leftmode = 1;
        else
            in_leftmode = 2;
        end
        
        in_ii=1; in_jj=1;
        while(in_ii <= intl_idxlen_2)
            if( X_new(in_ii)<=X_old(in_jj) )
                % less than left-reference
                if( X_new(in_ii)==X_new(in_ii+1))
                    in_ii = in_ii+1;
                else
                    if( X_old(in_jj)==X_old(in_jj+1) )
                        in_jj = in_jj+1;
                    else
                        if(X_new(in_ii)<X_old(in_jj))
                            in_leftmode = 2;
                        else
                            in_leftmode = 1;
                        end
                        
                        if(in_leftmode == 1)
                            if(in_jj>1)
                                in_jj = in_jj-1;  %decrease reference index on X_old
                            else
                                Y_new(in_ii)=Y_old(in_jj);
                                in_ii = in_ii+1;
                            end
                        else
                            Y_new(in_ii)=Y_old(in_jj);
                            in_ii = in_ii+1;  %increase reference index on X_new
                        end
                    end
                end
                
            else
                %within left-reference and right reference

                if( (X_new(in_ii)<=X_old(in_jj+1)) )
                    %weights
                    wgt_L = abs(1-(X_new(in_ii)-X_old(in_jj))/(X_old(in_jj+1)-X_old(in_jj)));
                    wgt_R = abs(1-(X_old(in_jj+1)-X_new(in_ii))/(X_old(in_jj+1)-X_old(in_jj)));
                    %linear interpolate
                    Y_new(in_ii) = Y_old(in_jj)*wgt_L + Y_old(in_jj+1)*wgt_R;

                    in_ii = in_ii+1;
                else
                    %greater than right reference due to exceeding ref points
                    if(in_jj < intl_idxlen_1-1)
                        in_jj = in_jj+1;    %increase reference index on X_old
                    else
                        Y_new(in_ii:end) = Y_old(end);
                        in_ii = intl_idxlen_2+1;
                    end
                    
                end
            end
        end
        
    end



end
%{
%==========================================================
    specu_nmd = specu - max(specu) + nmd_dB_offset;
    % Finding the edges of the CFBG spectrum
    % idx1 to idx2 is the extracted top flat region at db_thresh from max offsetted val
    idx_u1 = find(specu_nmd >= db_thresh,1,'first');
    idx_u2 = floor(length(specu_nmd)*3/5)-1+find(specu_nmd(floor(length(specu_nmd)*3/5):length(specu_nmd)) <= db_thresh,1,'first');
    idxlen_u = idx_u2-idx_u1 + 1;
    
    % creating common wlen array. Assumes wlen_u is monotonic-increasing
    itpl_upsample_factor = 4;
    wlen2_len = ceil(idxlen_u * itpl_upsample_factor);
    wlen_u2 = linspace(wlen_u(idx_u1), wlen_u(idx_u2), wlen2_len);
    wlen_u2_stepsize = (wlen_u(idx_u2)-wlen_u(idx_u1))/(wlen2_len-1);
    
    %Lin-interpolate specu to spec_u2 following wlen_u2
    spec_u2 = zeros(wlen2_len,1);
    ii=1; jj=idx_u1;
    while(ii <= wlen2_len)
        if( wlen_u2(ii)<wlen_u(jj) )
            % less than left-reference, somehow
            jj = jj-1;  %decrease reference index
        else
            %within left-reference and right reference
            if( (wlen_u2(ii)<=wlen_u(jj+1)) )
                %weights
                wgt_L = abs(1-(wlen_u2(ii)-wlen_u(jj))/(wlen_u(jj+1)-wlen_u(jj)));
                wgt_R = abs(1-(wlen_u(jj+1)-wlen_u2(ii))/(wlen_u(jj+1)-wlen_u(jj)));
                %linear interpolate
                spec_u2(ii) = specu_nmd(jj)*wgt_L + specu_nmd(jj+1)*wgt_R;
                
                ii = ii+1;
            else
                %greater than right reference due to exceeding ref points
                jj = jj+1;
            end
        end
    end
%}
%==========================================================