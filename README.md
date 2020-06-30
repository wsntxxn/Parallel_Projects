# Parallel programming projects

## Project 1

to compile and run programs in project 1, run this command: 
```
cd proj1
make
cd utils
make
cd ../
./run.sh
```
Before running `run.sh`, please make sure code, configuration and data files are on all hosts.

The default running processes is 4 for all tasks.
You can change this setting when runing single task with different processes. 

## Project 2

to compile and run programs in project 2, run this command: 
```
cd proj2
make
./run.sh
```

`run.sh` will run each project with 4 threads on single process (openmp version) as well as using 2 threads on 2 processes (openmp + mpi version).
You can change this setting by changing corresponding configuration in `run.sh`.

## Project 3

to compile and run program in project 3, simply run:
```
cd proj3
./run.sh
```
However, this program is based on Hadoop. Before running `run.sh`, you should make sure your hadoop is running normally (at least hdfs and yarn).

Besides, the default running script assumes that data is stored in `proj3/data` in plain text form.
You should first place all plain-text data (London2013.csv, Munbai2013.csv...) in `proj3/data`. 

Since the output is very long, I do not print to the screen but to a file (`./data_1-3/small_file_output` and `./data_1-3/big_file_output` respectively by default).