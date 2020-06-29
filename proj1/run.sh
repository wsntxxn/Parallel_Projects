echo -e "Runing 1-1"
mpirun -np 4 -f ./config ./1-1 100000

echo -e "Running 1-2"
if [ ! -d ./data_1-2 ]; then
    mkdir data_1-2
fi
echo -e "Multiplication of 2 1024 x 1024 matrices"
./utils/generate_matrix 1024 1024 0 0 data_1-2/data_1024 1
for i in {1..4}; do
    scp -q ./data_1-2/data_1024 slave$i:`pwd`/data_1-2/
done
mpirun -np 4 -f ./config ./1-2 ./data_1-2/data_1024 ./data_1-2/data_1024 ./data_1-2/matmul_result
echo -e "Convolution by im2col and matrix multiplication"
./utils/generate_matrix 1024 1024 4 4 ./data_1-2/matrix_row 0 ./data_1-2/matrix 0
./utils/generate_matrix 16 512 4 4 ./data_1-2/kernel_col_512 1 0 0
for i in {1..4}; do
    scp -q ./data_1-2/matrix_row ./data_1-2/kernel_col_512 slave$i:`pwd`/data_1-2/
done
mpirun -np 4 -f ./config ./1-2 ./data_1-2/matrix_row ./data_1-2/kernel_col_512 ./data_1-2/conv_result
# ./utils/resize_matrix ./data_1-2/conv_result 256 256 ./data_1-2/conv_result_resized
echo -e "Pooling by matrix multiplication"
./utils/generate_matrix 16 1 4 4 ./data_1-2/kernel_pooling 1 0 1
for i in {1..4}; do
    scp -q ./data_1-2/kernel_pooling slave$i:`pwd`/data_1-2/
done
mpirun -np 4 -f ./config ./1-2 ./data_1-2/matrix_row ./data_1-2/kernel_pooling ./data_1-2/pooling_result
./utils/resize_matrix ./data_1-2/pooling_result 256 256 ./data_1-2/pooling_result_resized

echo -e "Running 1-3"
echo -e "Word count of small files"
mpirun -np 4 -f ./config ./1-3-small ./data_1-3/small_file ./data_1-3/small_file_output
echo -e "Word count of big file"
mpirun -np 4 -f ./config ./1-3-big ./data_1-3/big_file/big_100.txt ./data_1-3/big_file_output