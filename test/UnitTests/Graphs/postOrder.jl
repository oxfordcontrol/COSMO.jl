workspace()


 function postOrder(par,child)
        N = length(par)
        order = zeros(Int64,N)
        root = find(x->x==0,par)[1]
        stack = Array{Int64,1}(0)
        iii = N
        push!(stack,root)
        while !isempty(stack)
            v = pop!(stack)
            order[v] = iii
            iii-=1
            push!(stack,child[v]...)
        end

        return order
    end


N = 6
par = [2;4;4;5;0;5]
child = [Array{Int64}(0) for i=1:N]
for iii=1:N
    par[iii] != 0 && (push!(child[par[iii]],iii))
end


postOrder(par,child)