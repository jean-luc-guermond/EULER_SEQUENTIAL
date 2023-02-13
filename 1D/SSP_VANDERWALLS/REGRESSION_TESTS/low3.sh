run="srun-serial ../EXECUTABLE/ssp_euler_1D_vanderwalls.exe regression"

cp data_low_order_case3 data
cp mesh_case3.FEM mesh.FEM
cp reference_regression_test_low_order_case3 reference_regression_test
$run
