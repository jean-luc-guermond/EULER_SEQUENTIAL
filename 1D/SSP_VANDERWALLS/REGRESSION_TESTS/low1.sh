run="srun-serial ../EXECUTABLE/ssp_euler_1D_vanderwalls.exe regression"

cp data_low_order_case1 data
cp mesh_case1.FEM mesh.FEM
cp reference_regression_test_low_order_case1 reference_regression_test
$run
