git clone https://github.com/oxfordcontrol/osqp_benchmarks/
cd osqp_benchmarks
git checkout b62d31a
cd ..
mv osqp_benchmarks/problem_classes/maros_meszaros_data/*.mat ./
rm -rf osqp_benchmarks