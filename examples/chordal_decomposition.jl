using COSMO, LinearAlgebra, JuMP, Random, SparseArrays

# Define problem matrices
n = 9
m = 2

A1 = [-4.0 0.0 -2.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 -3.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; -2.0 -1.0 -2.0 0.0 0.0 5.0 4.0 -4.0 0.0; 0.0 0.0 0.0 -4.0 -5.0 0.0 0.0 3.0 0.0; 0.0 0.0 0.0 -5.0 4.0 0.0 0.0 2.0 0.0; -1.0 0.0 5.0 0.0 0.0 5.0 -4.0 -4.0 -5.0; 0.0 0.0 4.0 0.0 0.0 -4.0 -1.0 -1.0 -3.0; 0.0 0.0 -4.0 3.0 2.0 -4.0 -1.0 2.0 -2.0; 0.0 0.0 0.0 0.0 0.0 -5.0 -3.0 -2.0 -3.0];
A2 = [-5.0 0.0 3.0 0.0 0.0 -2.0 0.0 0.0 0.0; 0.0 -3.0 -5.0 0.0 0.0 0.0 0.0 0.0 0.0; 3.0 -5.0 3.0 0.0 0.0 5.0 -4.0 -5.0 0.0; 0.0 0.0 0.0 3.0 2.0 0.0 0.0 -2.0 0.0; 0.0 0.0 0.0 2.0 4.0 0.0 0.0 -3.0 0.0; -2.0 0.0 5.0 0.0 0.0 1.0 -5.0 -2.0 -4.0; 0.0 0.0 -4.0 0.0 0.0 -5.0 -2.0 -3.0 3.0; 0.0 0.0 -5.0 -2.0 -3.0 -2.0 -3.0 5.0 3.0; 0.0 0.0 0.0 0.0 0.0 -4.0 3.0 3.0 -4.0];
B= [-0.11477375644968069 0.0 6.739182490600791 0.0 0.0 -1.2185593245043502 0.0 0.0 0.0; 0.0 1.2827680528587497 -5.136452036888789 0.0 0.0 0.0 0.0 0.0 0.0; 6.739182490600791 -5.136452036888789 7.344770673489607 0.0 0.0 -0.2224400187044442 -10.505300166831221 -1.2627361794562273 0.0; 0.0 0.0 0.0 10.327710040060499 8.91534585379813 0.0 0.0 -6.525873789637007 0.0; 0.0 0.0 0.0 8.91534585379813 0.8370459338528677 0.0 0.0 -6.210900615408826 0.0; -1.2185593245043502 0.0 -0.2224400187044442 0.0 0.0 -3.8185953011245024 -0.994033914192722 2.8156077981712997 1.4524716674219218; 0.0 0.0 -10.505300166831221 0.0 0.0 -0.994033914192722 0.029162208619863517 -2.8123790276830745 7.663416446183705; 0.0 0.0 -1.2627361794562273 -6.525873789637007 -6.210900615408826 2.8156077981712997 -2.8123790276830745 4.71893305728242 6.322431630550857; 0.0 0.0 0.0 0.0 0.0 1.4524716674219218 7.663416446183705 6.322431630550857 0.5026094532322212];
c = [-0.21052661285686525, -1.263324575834677];

# Solve with JuMP + COSMO


# Solve without clique merging
model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer, "merge_strategy" => COSMO.NoMerge, "complete_dual" => true));
@variable(model, x[1:m]);
@objective(model, Min, c' * x )
@constraint(model, Symmetric(B - A1  .* x[1] - A2 .* x[2] )  in JuMP.PSDCone());
JuMP.optimize!(model)


# Solve with a clique graph based strategy
model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer, "merge_strategy" => COSMO.CliqueGraphMerge, "complete_dual" => true));
@variable(model, x[1:m]);
@objective(model, Min, c' * x )
@constraint(model, Symmetric(B - A1  .* x[1] - A2 .* x[2] )  in JuMP.PSDCone());
JuMP.optimize!(model)


# Solve with a parent-child merge strategy
pc_strategy = with_options(COSMO.ParentChildMerge, t_size = 3, t_fill = 3)
model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer, "merge_strategy" => pc_strategy, "complete_dual" => true));
@variable(model, x[1:m]);
@objective(model, Min, c' * x )
@constraint(model, Symmetric(B - A1  .* x[1] - A2 .* x[2] )  in JuMP.PSDCone());
JuMP.optimize!(model)
