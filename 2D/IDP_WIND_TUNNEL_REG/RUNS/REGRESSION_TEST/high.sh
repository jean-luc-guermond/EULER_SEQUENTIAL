cp data_high_order data
cp reference_regression_test_high_order reference_regression_test
srun-serial ../../EXECUTABLE/idp_euler.exe regression
rm reference_regression_test
