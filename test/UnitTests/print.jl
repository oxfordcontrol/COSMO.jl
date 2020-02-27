using COSMO, Test, Random

ws = COSMO.Workspace()
iter = 10
cost = 20.
r_prim = 1.5e-3
r_dual = 1.2e-2
status = :Solved
rt = 0.7
time = 0.1


@testset "Printing" begin
   @test COSMO.print_iteration(ws, iter, cost, r_prim, r_dual, time) == nothing
   @test COSMO.print_iteration(ws, ws.settings.check_termination, cost, r_prim, r_dual, time) == nothing
   @test COSMO.print_result(status, iter, cost, rt) == nothing
end
nothing
