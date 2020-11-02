function [output] = fn_least_squares(x_in, y_in, order, output_mode)
    % Author: Alex Wong (almw994@uowmail.edu.au)
    % fn_least_squares v1.0 25-07-2019
    %
    % fn_least_squares calculates the Nth-degree least square polynomial
    % of an x-y data series for best-fitting
    % syntax:
    % fn_least_squares(x_in, y_in, order, output_mode)
    % output_mode: 
    %   mode 0 - coefficients only 
    %       Eg. order 1 (linear) - output = [gradient,intercept]
    %   mode 1 - polynomial vector output 
    
    if( max(size(x_in) ~= size(y_in)) || ( min(size(x_in))~= 1 ))
        error('x and y inputs are either not vectors or agreeing in dimensions');
    end
    
    dat_len = length(x_in);
    N = order;
    
    if(mod(N,1)~=0 || N<1)
        error('polynomial order must be integer >= 1');
    end
    
    if(output_mode ~=0 && output_mode ~=1)
        error('output_mode must be 0 or 1. See syntax');
    end
    
    if(dat_len <= N)
        error('input length must be greater than order');
    end

    
    fit_A = zeros(N+1,N+1);
    %temp_A = fit_A;
    fit_Y = zeros(N+1,1);
    fit_coeffs = zeros(N+1,1);
    fit_result = zeros(dat_len,1);
        
    %creating S_1_0 to S_M_0
    fit_A(1,1)=dat_len;
    ii=2; jj=1; kk=1;
    while(ii<=N+1)
        kk=1;
        while(kk<=dat_len)
            fit_A(ii,jj) = fit_A(ii,jj) + x_in(kk)^(ii+jj-2);
            kk=kk+1;
        end
        ii=ii+1;
    end
    ii=ii-1; jj=jj+1;
    while(jj<=N+1)
        kk=1;
        while(kk<=dat_len)
            fit_A(ii,jj) = fit_A(ii,jj) + x_in(kk)^(ii+jj-2);
            kk=kk+1;
        end        
        jj=jj+1;
    end
    
    %filling fit_A array of S_M_0
    ii=2; jj=1;
    while(ii<=N+1)
        jj=1; kk=ii;
        while(kk>=2)
            fit_A(kk-1,jj+1) = fit_A(kk,jj);
            
            jj=jj+1;
            kk=kk-1;
        end
        ii = ii+1;
    end
    ii=ii-1; jj=2;
    while(jj<=N+1)
        kk=ii; hh=jj; 
        while(hh<=N)
            fit_A(kk-1,hh+1) = fit_A(kk,hh);
            
            hh=hh+1;
            kk=kk-1;
        end
        jj = jj+1;
    end
   
    %fit_Y vector 
    ii=1;
    while(ii<=N+1)
        kk=1;
        while(kk<=dat_len)
            fit_Y(ii) = fit_Y(ii) + x_in(kk)^(ii-1)*y_in(kk);
            kk = kk+1;
        end
        ii=ii+1;
    end
    
    
    %fit_A = fit_A'
    
    %diagonally mirroring (top left exchanged with bottom right) of fit_A
    %{
    ii=1; jj=1;
    temp_A = fit_A;
    while(ii<= N+1)
        jj=1;
        while(jj<=N+1)
            fit_A(ii,jj)=temp_A(N+2-ii,N+2-jj);            
            jj = jj+1;
        end
        ii = ii+1;
    end
    %}
    fit_A = fliplr(flipud(fit_A))';
    fit_Y = flipud(fit_Y);
    
    fit_coeffs = linsolve(fit_A,fit_Y);
    
    fit_coeffs = flipud(fit_coeffs);
    fit_coeffs = fit_coeffs';
    
    ii=1;
    while(ii<=dat_len)
        jj=1;
        while(jj<=N+1)
            fit_result(ii) = fit_result(ii) + fit_coeffs(jj)*x_in(ii)^(jj-1);
            jj = jj+1;
        end
        ii=ii+1;
    end
    
           
    switch output_mode
        case 0
            output = fit_coeffs;
        case 1
            output = fit_result;
        otherwise
            error('error in switch expression');
    end
    
    
    %{
    cfit_A = [...
            [S120 S110 S10 S90 S80 S70 S60]; ...
            [S110 S100 S90 S80 S70 S60 S50]; ...
            [S100 S90 S80 S70 S60 S50 S40]; ...
            [S90 S80 S70 S60 S50 S40 S30]; ...
            [S80 S70 S60 S50 S40 S30 S20]; ...
            [S70 S60 S50 S40 S30 S20 S10]; ...
            [S60 S50 S40 S30 S20 S10 S00] ...
            ];
    cfit_Y = [S61; S51; S41; S31; S21; S11; S01];
    %}
    
%{    

if(do_lsq==1)
    x_cfit = -5:0.02:xlen+5;
    ii = 1;
    S00 = x_idxlen;
    S10 = 0;
    S20 = 0;
    S30 = 0;
    S40 = 0;
    S50 = 0;
    S60 = 0;
    S70 = 0;
    S80 = 0;
    S90 = 0;
    S100 = 0;
    S110 = 0;
    S120 = 0;
    S01 = 0;
    S11 = 0;
    S21 = 0;
    S31 = 0;
    S41 = 0;
    S51 = 0;
    S61 = 0;

    while(ii <= x_idxlen)
        S10 = S10 + x_sLR_ST(ii);
        S20 = S20 + x_sLR_ST(ii)^2;
        S30 = S30 + x_sLR_ST(ii)^3;
        S40 = S40 + x_sLR_ST(ii)^4;
        S50 = S50 + x_sLR_ST(ii)^5;
        S60 = S60 + x_sLR_ST(ii)^6;
        S70 = S70 + x_sLR_ST(ii)^7;
        S80 = S80 + x_sLR_ST(ii)^8;
        S90 = S90 + x_sLR_ST(ii)^9;
        S100 = S100 + x_sLR_ST(ii)^10;
        S110 = S110 + x_sLR_ST(ii)^11;
        S120 = S120 + x_sLR_ST(ii)^12;
        S01 = S01 + (x_sLR_ST(ii)^0)*combined_3sLR_ST_2(ii);
        S11 = S11 + (x_sLR_ST(ii)^1)*combined_3sLR_ST_2(ii);
        S21 = S21 + (x_sLR_ST(ii)^2)*combined_3sLR_ST_2(ii);
        S31 = S31 + (x_sLR_ST(ii)^3)*combined_3sLR_ST_2(ii);
        S41 = S41 + (x_sLR_ST(ii)^4)*combined_3sLR_ST_2(ii);
        S51 = S51 + (x_sLR_ST(ii)^5)*combined_3sLR_ST_2(ii);
        S61 = S61 + (x_sLR_ST(ii)^6)*combined_3sLR_ST_2(ii);

        ii = ii+1;
    end

    switch(lsq_order)
        case {2}
            cfit_A = [...
                    [S40 S30 S20]; ...
                    [S30 S20 S10]; ...
                    [S20 S10 S00]
                    ];
            cfit_Y = [S21; S11; S01];
            cfit_coeffs = linsolve(cfit_A,cfit_Y);

            y_cfit = cfit_coeffs(1).*x_cfit.^2 ...
                + cfit_coeffs(2).*x_cfit + cfit_coeffs(3);

        case {3}
            cfit_A = [...
                    [S60 S50 S40 S30]; ...
                    [S50 S40 S30 S20]; ...
                    [S40 S30 S20 S10]; ...
                    [S30 S20 S10 S00]
                    ];
            cfit_Y = [S31; S21; S11; S01];
            cfit_coeffs = linsolve(cfit_A,cfit_Y);

            y_cfit = cfit_coeffs(1).*x_cfit.^3 + cfit_coeffs(2).*x_cfit.^2 ...
                + cfit_coeffs(3).*x_cfit + cfit_coeffs(4);

        case {4}
            cfit_A = [...
                    [S80 S70 S60 S50 S40]; ...
                    [S70 S60 S50 S40 S30]; ...
                    [S60 S50 S40 S30 S20]; ...
                    [S50 S40 S30 S20 S10]; ...
                    [S40 S30 S20 S10 S00] ...
                    ];
            cfit_Y = [S41; S31; S21; S11; S01];
            cfit_coeffs = linsolve(cfit_A,cfit_Y);

            y_cfit = cfit_coeffs(1).*x_cfit.^4 + cfit_coeffs(2).*x_cfit.^3 + cfit_coeffs(3).*x_cfit.^2 ...
                + cfit_coeffs(4).*x_cfit + cfit_coeffs(5);


        case {5}
            cfit_A = [...
                    [S100 S90 S80 S70 S60 S50]; ...
                    [S90 S80 S70 S60 S50 S40]; ...
                    [S80 S70 S60 S50 S40 S30]; ...
                    [S70 S60 S50 S40 S30 S20]; ...
                    [S60 S50 S40 S30 S20 S10]; ...
                    [S50 S40 S30 S20 S10 S00] ...
                    ];
            cfit_Y = [S51; S41; S31; S21; S11; S01];
            cfit_coeffs = linsolve(cfit_A,cfit_Y);

            y_cfit = cfit_coeffs(1).*x_cfit.^5 + cfit_coeffs(2).*x_cfit.^4 + cfit_coeffs(3).*x_cfit.^3 + cfit_coeffs(4).*x_cfit.^2 ...
                + cfit_coeffs(5).*x_cfit + cfit_coeffs(6);


        case {6}
            cfit_A = [...
                    [S120 S110 S10 S90 S80 S70 S60]; ...
                    [S110 S100 S90 S80 S70 S60 S50]; ...
                    [S100 S90 S80 S70 S60 S50 S40]; ...
                    [S90 S80 S70 S60 S50 S40 S30]; ...
                    [S80 S70 S60 S50 S40 S30 S20]; ...
                    [S70 S60 S50 S40 S30 S20 S10]; ...
                    [S60 S50 S40 S30 S20 S10 S00] ...
                    ];
            cfit_Y = [S61; S51; S41; S31; S21; S11; S01];
            cfit_coeffs = linsolve(cfit_A,cfit_Y);

            y_cfit = cfit_coeffs(1).*x_cfit.^6 + cfit_coeffs(2).*x_cfit.^5 + cfit_coeffs(3).*x_cfit.^4 + cfit_coeffs(4).*x_cfit.^3 + cfit_coeffs(5).*x_cfit.^2 ...
                + cfit_coeffs(6).*x_cfit + cfit_coeffs(7);



    end
%}
    
    
    
    
    
end

