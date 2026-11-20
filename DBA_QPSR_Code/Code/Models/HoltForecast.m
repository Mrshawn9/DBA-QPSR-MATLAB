function y = HoltForecast(x, A, H)

A = A(:)';
n = numel(A);

para = x(1);
para1 = x(2);
para2 = x(3);

s1 = zeros(1, n);
s2 = zeros(1, n);
s3 = zeros(1, n);

for k = 1:n
    if k == 1
        s1(k) = A(k);
    else
        vals = zeros(1, k);
        vals(1) = ((1-para)^(k-1))*A(1);

        for j = 2:k
            vals(j) = para*((1-para)^(k-j))*A(j);
        end

        s1(k) = sum(vals);
    end
end

for k = 1:n
    if k == 1
        s2(k) = s1(k);
    else
        vals = zeros(1, k);
        vals(1) = ((1-para1)^(k-1))*s1(1);

        for j = 2:k
            vals(j) = para1*((1-para1)^(k-j))*s1(j);
        end

        s2(k) = sum(vals);
    end
end

for k = 1:n
    if k == 1
        s3(k) = s2(k);
    else
        vals = zeros(1, k);
        vals(1) = ((1-para2)^(k-1))*s2(1);

        for j = 2:k
            vals(j) = para2*((1-para1)^(k-j))*s2(j);
        end

        s3(k) = sum(vals);
    end
end

a = 3*s1-3*s2+s3;
b = (para/2*((1-para)^2))* ...
    ((6-5*para)*s1-2*(5-4*para)*s2+(4-3*para)*s3);
c = ((para/2*(1-para))^2)*(s1-2*s2+s3);

y = nan(1, H);

for m = 1:H
    y(m) = a(end)+m*b(end)+(m^2)*c(end);
end

if any(~isfinite(y)) || ~isreal(y)
    y(:) = NaN;
end
end
