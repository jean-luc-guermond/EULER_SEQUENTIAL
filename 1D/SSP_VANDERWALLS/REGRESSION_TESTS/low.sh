run="srun-serial ../../EXECUTABLE/ssp_euler_1D_vanderwalls.exe regression"

cp data_low_order_case0 data
cp mesh_case0.FEM mesh.FEM
cp reference_regression_test_low_order_case0 reference_regression_test
$run

#sleep 1

#cp data_low_order_case1 data
#cp mesh_case1.FEM mesh.FEM
#cp reference_regression_test_low_order_case1 reference_regression_test
#$run

#sleep 1

#cp data_low_order_case2 data
#cp mesh_case2.FEM mesh.FEM
#cp reference_regression_test_low_order_case2 reference_regression_test
#$run

#sleep 1

#cp data_low_order_case3 data
#cp mesh_case3.FEM mesh.FEM
#cp reference_regression_test_low_order_case3 reference_regression_test
#$run
