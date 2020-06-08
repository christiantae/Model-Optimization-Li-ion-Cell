function [x_best, curr_cost] = grad_descent(x0, dataset)
learning_rate = 0.1;
num_iterations = 10;
epsilon = 1e-8
x = x0'; 
x_best = 0;
cost = Inf;
for n = 1 : num_iterations
    grads = estimate_derivative(x, dataset);
    %d = -grads/(norm(x) + epsilon)
    d = -grads 
    x =  x + learning_rate.*d;
    curr_cost = ECM_fit(x, dataset);
    if (curr_cost < cost)
        x_best = x;
        cost = curr_cost;
    end
end

end 