function f = HoltObjective(x, A)

A = A(:)';
n = numel(A);

q = zeros(1, n);
p = zeros(1, n);
qq = zeros(1, n);

for k = 1:n
    if k == 1
        q(k) = A(k);
    else
        vals = zeros(1, k);
        vals(1) = ((1-x(1))^(k-1))*A(1);

        for j = 2:k
            vals(j) = x(1)*((1-x(1))^(k-j))*A(j);
        end

        q(k) = sum(vals);
    end
end

for k = 1:n
    if k == 1
        p(k) = q(k);
    else
        vals = zeros(1, k);
        vals(1) = ((1-x(2))^(k-1))*q(1);

        for j = 2:k
            vals(j) = x(2)*((1-x(2))^(k-j))*q(j);
        end

        p(k) = sum(vals);
    end
end

for k = 1:n
    if k == 1
        qq(k) = p(k);
    else
        vals = zeros(1, k);
        vals(1) = ((1-x(3))^(k-1))*p(1);

        for j = 2:k
            vals(j) = x(3)*((1-x(3))^(k-j))*p(j);
        end

        qq(k) = sum(vals);
    end
end

x1 = 3*q-3*p+qq;
f = sum(abs((x1-A)./A));

if ~isfinite(f) || ~isreal(f)
    f = 1e30;
end
end
