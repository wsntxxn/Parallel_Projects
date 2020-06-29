#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <unistd.h>
#include <string.h>
#include <mpi.h>

#include <unordered_map>
#include <iostream>

#define EMPTY_MSG 0
#define FILE_MSG 1
#define RESULT_MSG 2

#define MAX_LINE_LENGTH 1024
#define MAX_NUMBER_LENGTH 5
#define MAX_STR_LENGTH 1048000

#define DELIMITERS " \t\n\r\f\v"

using namespace std;

void serialize(unordered_map<string, int> map, char** result) {
    char *str = (char *) malloc(MAX_STR_LENGTH * sizeof(char));
    str[0] = '\0';

	unordered_map<string, int>::iterator iter;
    char *count = (char *) malloc(MAX_NUMBER_LENGTH * sizeof(char));
	for(iter = map.begin(); iter != map.end(); iter++) {
        strcat(str, ((iter -> first).c_str()));
        sprintf(count, "%d", iter -> second);
        strcat(str, " ");
        strcat(str, count);
        strcat(str, " ");
    }

    strcpy(*result, str);
    free(str);
    free(count);
}

void deserialize(char *result, unordered_map<string, int> *map) {

    char *pch;
    pch = strtok (result, " ");

    int idx = 0;
    char *word;
    int count;

    word = (char *) malloc(MAX_STR_LENGTH * sizeof(char));
    while (pch != NULL)
    {
        if (idx % 2 == 0) // word
            strcpy(word, pch);
        else {// count
            count = atoi(pch);
            (*map)[word] = count;
        }

        pch = strtok (NULL, " ");
        idx++;
    }

    free(word);

}

unordered_map<string, int> build_wcmap_from_file(char* filename) {

    unordered_map<string, int> map; 
    FILE *fp;
    int line_length;
    char line[MAX_LINE_LENGTH];
    
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Opening %s error\n", filename);
        return map;
    }

    while(fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        line_length = strlen(line);
        line[line_length - 1] = '\0';
        line_length--;
        char *word;
        word = strtok(line, DELIMITERS);
        while(word != NULL)
        {
            if (map.count(word) == 0)
                map[word] = 1;
            else
                map[word] += 1;
            //printf("%s\t", word);
            word = strtok(NULL, DELIMITERS);
        }
        //printf("\n");
    }

    return map;
    
}


void manager(char* argv[], int processes) {
    
    DIR *dp;
    struct dirent *dirp;
    char cwd[80];
    int terminated;
    
    getcwd(cwd, sizeof(cwd));
    chdir(argv[1]);
    dp = opendir("./");
    if(dp == NULL) {
        printf("Opening %s error\n", argv[1]);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    terminated = 0;
    char *result_part;
    unordered_map<string, int> map_part;
    unordered_map<string, int> map_global;
    int total_count = 0;
    unordered_map<string, int>::iterator iter;
    MPI_Status status;
    int src;
    int tag;
    char *filename;

    unordered_map<string, int> rank2filename;

    result_part = (char *) malloc(MAX_STR_LENGTH * sizeof(char));

    int file_idx = 0;

    do {            
        MPI_Recv(result_part, MAX_STR_LENGTH, MPI_CHAR, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        src = status.MPI_SOURCE;
        tag = status.MPI_TAG;

        int length;
        MPI_Get_count(&status, MPI_CHAR, &length);
        result_part[length] = '\0';
        
        dirp = readdir(dp);

        if (tag == RESULT_MSG) {

            if (dirp != NULL) { // still have remaining files 
                filename = dirp -> d_name;
                MPI_Send(filename, strlen(filename), MPI_CHAR, src, FILE_MSG, MPI_COMM_WORLD);
            }
            else { // this process finishes counting
                MPI_Send (NULL, 0, MPI_INT, src, EMPTY_MSG, MPI_COMM_WORLD);
                terminated++;
            }
            
            // deserialize result from this file and update global result
            deserialize(result_part, &map_part);


            for(iter = map_part.begin(); iter != map_part.end(); iter++) {
                if (map_global.count(iter -> first) == 0)
                    map_global[iter -> first] = iter -> second;
                else
                    map_global[iter -> first] += iter -> second;
                total_count += iter -> second;

            }

            map_part.clear();

        }
        else { // EMPTY MSG, means initial request for work
            filename = dirp -> d_name;            
            MPI_Send(filename, strlen(filename), MPI_CHAR, src, FILE_MSG, MPI_COMM_WORLD);
        }

    } while(terminated < processes - 1);
    
    char *word;
    FILE *fout;
    chdir(cwd);
    word = (char *) malloc(sizeof(char) * MAX_STR_LENGTH);
    fout = fopen(argv[2], "w");
    fprintf(fout, "word\t--\tfrequency\n");
    for(iter = map_global.begin(); iter != map_global.end(); iter++) {
        strcpy(word, (iter -> first).c_str());
        fprintf(fout, "%s\t--\t%f\n", word, (float)(iter -> second) / total_count);
    }
    fclose(fout);

}

void worker(char *argv[], int processes, int rank) {
    chdir(argv[1]);

    MPI_Request pending;
    MPI_Status status;

    char *result;
    char *filename;
    int filename_length;
    unordered_map<string, int> map;


    result = (char *) malloc(sizeof(char) * MAX_STR_LENGTH);
    // Initial request for work
    MPI_Isend(NULL, 0, MPI_CHAR, 0, EMPTY_MSG, MPI_COMM_WORLD, &pending);

    while (1) {
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &filename_length);

        filename = (char *) malloc(sizeof(char) * filename_length);

        if (status.MPI_TAG == FILE_MSG) {
            MPI_Recv(filename, filename_length, MPI_CHAR, 0, FILE_MSG, MPI_COMM_WORLD, &status);
            filename[filename_length] = '\0';

            map = build_wcmap_from_file(filename);
            
            unordered_map<string, int>::iterator iter;

            serialize(map, &result);
            MPI_Send(result, strlen(result), MPI_CHAR, 0, RESULT_MSG, MPI_COMM_WORLD);

        }
        else { // tag == EMPTY MSG
            break;
        }
    }
}


int main(int argc, char* argv[]) {
    int rank;
    int processes;

    MPI_Comm worker_comm;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    if (argc != 3) {
        if (!rank) {
            printf("Usage: %s <input directory> <output filename>\n", argv[0]);

        }
    }
    else if (processes < 2) {
        printf("This program needs at least two processes\n");
    }
    else {
        MPI_Barrier(MPI_COMM_WORLD);
        double elapsed_time = -MPI_Wtime();

        if (!rank) {
            manager(argv, processes);
        }
        else {
            worker(argv, processes, rank);
        }

        elapsed_time += MPI_Wtime();
        if (!rank) {
            printf("elapsed time: %lfs\n", elapsed_time);
        }
    }

    MPI_Finalize();
    return 0;
}
