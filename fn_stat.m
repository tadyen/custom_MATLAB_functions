function output = fn_stat(data_in, method, order)
    % Author: Alex Wong (almw994@uowmail.edu.au)
    % fn_stat v1.0 25-07-2019
    %
    % calculates mean, variance and std_deviation of a sample
    % syntax: fn_stat(dat_in, mode, order)
    % methods:
    %   'p' - Population. Normal calculation of mean, variance, and stddev
    %         of a 1D array of values.
    %         output = [mean, variance, std_dev];
    %
    %   's' - Sample. Similar to 'p' but calculation is done using 'N-1' 
    %         instead of 'N'.
    %         output = [mean, variance, std_dev];
    %
    %   'r' - Regression. Uses fn_least_squares for regression.
    %       Takes an (n x 2) vector for x and y values, and calculates 
    %       against a polynomial regression instead of the sample mean. 
    %       The 'order' argument is used for the order of the regression 
    %       polynomial. Sample (N-1) method
    %       is used. Also detects outliers by finding the relative data
    %       index of the largest contributor to the variance. Output is a
    %       struct with the following:
    %        ===
    %       output.mean (of sample)
    %       output.variance (of sample vs regression values)
    %       output.stddev (of sample vs regression values)
    %       output.vects.variance   (variances of each data)
    %       output.vects.stddev     (stddev of each data)
    %       output.vects.score      (stddevs normalised to stddev)
    %       output.vects.regression (regression line)
    %       output.outlier_index      (relative index of outlier)
    %       output.outlier_score    (specific_stddev proportion to stddev)
    %        ===
    
    switch method
        case 'p'
            ;
        case 's'
            ;
        case 'r'
            ;
        otherwise
            error('Method argument is incorrect. See syntax');
    end
    
    data_size = size(data_in);
    if(method == 'r')
        if(data_size(2)~=2)
            error('fn_stat input is incorrect for r method');
        end
        data_len = data_size(1);
    else
        if(min(size(data_in))~=1)
            error('fn_stat input not a 1D array');
        end    

        data_len = length(data_in);
        if(data_len == 1)
            error('fn_stat input not a vector');
        end
    end
    

    mean = 0;
    variance = 0;
    stddev = 0;
    
    ii=1;
    while(ii <= data_len)
        if(method ~= 'r')
            mean = mean + data_in(ii); 
        else
            mean = mean + data_in(ii,2);
        end
        ii = ii+1;
    end
    mean = mean/data_len;
    
    %calculate regression
    if(method=='r')
        xx = data_in(:,1);
        yy = data_in(:,2);
        zz = fn_least_squares(xx,yy,order,1);
    end
    
    
    
    switch method
        case 'p'
            ii=1;
            while(ii <= data_len)
                variance = variance + (data_in(ii) - mean)^2; 
                ii = ii+1;
            end
            variance = variance/(data_len);
            stddev = sqrt(variance);
        case 's'
            ii=1;
            while(ii <= data_len)
                variance = variance + (data_in(ii) - mean)^2; 
                ii = ii+1;
            end
            variance = variance/(data_len-1);
            stddev = sqrt(variance);
        case 'r'
            output.vects.variance=zeros(data_len,1);
            output.vects.stddev=zeros(data_len,1);
            ii=1;
            while(ii <= data_len)
                output.vects.variance(ii) = (yy(ii) - zz(ii))^2/(data_len-1);
                output.vects.stddev(ii) = sqrt(output.vects.variance(ii));
                variance = variance + (yy(ii) - zz(ii))^2/(data_len-1); 
                ii = ii+1;
            end
            stddev = sqrt(variance);
            output.mean = mean;
            output.variance = variance;
            output.stddev = stddev;
            output.vects.regression = zz;
            output.vects.score = output.vects.stddev/stddev;
            output.outlier_index = find(output.vects.score==max(output.vects.score));
            output.outlier_score = output.vects.score(output.outlier_index);
        otherwise
            error('Method argument is incorrect. See syntax');
    end
    
    if(isreal(mean)==false || isreal(variance)==false || isreal(stddev)==false)
        error('complex output found');
    end
    
    if(method ~= 'r')
        output = [mean, variance, stddev];
    else
        ;
    end
    
end

