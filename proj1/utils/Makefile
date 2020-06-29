ALL=generate_matrix display_matrix resize_matrix

all: ${ALL}

generate_matrix: generate_matrix.c
	gcc -o $@ generate_matrix.c

display_matrix: display_matrix.c
	gcc -o $@ display_matrix.c

resize_matrix: resize_matrix.c
	gcc -o $@ resize_matrix.c

clean:
	rm -f ${ALL}