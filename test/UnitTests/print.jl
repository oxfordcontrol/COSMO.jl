using COSMO, Test, Random

# This test is precision agnostic and will be run with T precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end

@testset "Printing" begin
for T in UnitTestFloats

   ws = COSMO.Workspace{T}()
   iter = 10
   cost = 20.
   r_prim = 1.5e-3
   r_dual = 1.2e-2
   status = :Solved
   rt = 0.7
   safeguard = true
   safeguarding_iter = 1

   @testset "Printing (T = $(T))" begin
      @test COSMO.print_iteration(ws, iter, cost, r_prim, r_dual) == nothing
      @test COSMO.print_iteration(ws, ws.settings.check_termination, cost, r_prim, r_dual) == nothing
      @test COSMO.print_result(status, iter, safeguarding_iter, cost, rt, safeguard) == nothing
   end
end
end
nothing
