function X = initialization(SearchAgents_no, dim, ub, lb)

Boundary_no = size(ub, 2);

if Boundary_no == 1
    X = rand(SearchAgents_no, dim).*(ub-lb)+lb;
end

if Boundary_no > 1
    for i = 1:dim
        ub_i = ub(i);
        lb_i = lb(i);
        X(:, i) = rand(SearchAgents_no, 1).*(ub_i-lb_i)+lb_i;
    end
end
