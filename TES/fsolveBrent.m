function [s,f_s,count] = fsolveBrent(f,a,b,tol_y)
    %FSOLVEBRENT Using Brent's method solve f(x) = 0;
    % must have f(a)*f(b) < 0
    %   Detailed explanation goes here
    max_count = 1000;
    count = 0;
    if nargin == 3
        tol_y = 1e-6;
    end
    f_a = f(a);
    if abs(f_a) < tol_y
        s = a; f_s = f_a; 
        return
    end
    f_b = f(b);
    if abs(f_b) < tol_y
        s = b; f_s = f_b;
        return
    end
    if f_a*f_b > 0
        error('Root is not bracketed')
    end
    tol_x = max(max(a/abs(f_a),b/abs(f_b))*tol_y,1e-9); 
    if abs(f_a) < abs(f_b)
        temp = a;
        a = b;
        b = temp;
        temp = f_a;
        f_a = f_b;
        f_b = temp;
    end
    c = a;
    f_c = f_a;
    mflag = true;
    while abs(f_b) > tol_y && count < max_count
        count = count+1;
        if f_a ~= f_c && f_b ~= f_c
            s = a*f_b*f_c/((f_a-f_b)*(f_a-f_c))+b*f_a*f_c/((f_b-f_a)*(f_b-f_c))+c*f_a*f_b/((f_c-f_a)*(f_c-f_b));
        else
            s = b-f_b*(b-a)/(f_b-f_a);
        end
        bounds = sort([(3*a+b)/4,b]);
        if ( s <= bounds(1) || s >= bounds(2)) || ...
                    ( mflag && ((abs(s-b) >= abs(b-c)/2) || (abs(b-c) < tol_x))) || ...
                    (~mflag && ((abs(s-b) >= abs(c-d)/2) || (abs(c-d) < tol_x)))

            s = (a+b)/2;
            mflag = true;
        else
            mflag = false;
        end
        f_s = f(s);
        d = c;
        c = b;
        f_c = f_b;
        if f_a*f_s < 0
            b = s;
            f_b = f_s;
        else
            a = s;
            f_a = f_s;
        end
        if abs(f_a) < abs(f_b)
            temp = a;
            a = b;
            b = temp;
            temp = f_a;
            f_a = f_b;
            f_b = temp;
        end
    end

end

