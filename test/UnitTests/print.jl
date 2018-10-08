using COSMO, Test, Random

settings = COSMO.Settings()
iter = 10
cost = 20.
r_prim = 1.5e-3
r_dual = 1.2e-2
status = :Solved
rt = 0.7


@testset "Printing" begin
   @test COSMO.Printing.printIteration(settings,iter,cost,r_prim,r_dual) == nothing
   @test COSMO.Printing.printIteration(settings,settings.check_termination,cost,r_prim,r_dual) == nothing
   @test COSMO.Printing.printResult(status,iter,cost,rt) == nothing
end
nothing
