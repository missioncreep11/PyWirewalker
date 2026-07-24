function Coeffs = FitKolmogorov(f,y)
    Afunc = @(fin,yin) [fin./fin fin.^(-5/3)];

    if size(f,1)==1
        f = f';
    end
    if size(y,1)==1
        y = y';
    end
    
    mask = (isnan(f) | isnan(y));
    f(mask) = [];
    y(mask) = [];
    
    A = Afunc(f,y);

    Coeffs = (A'*A)^-1*A'*y;
end