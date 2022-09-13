#ifndef INV_INPUT_H
#define INV_INPUT_H


/**
 * @brief A struct holding a sparse CSR matrix and its adjacency structure
 * 
 */
typedef struct InputCSR {
    int n;
    int m;
    int* i;
    int* j;
    double* val;
    int* nnz;
    int* xadj;
    int* yadj;
} InputCSR;


/**
 * @brief A constructor for InputCSR
 * 
 * @param n 
 * @param m 
 * @return InputCSR* 
 */
InputCSR* input_create(int n, int m);


/**
 * @brief A destructor for InputCSR
 * 
 * @param input 
 */
void input_destroy(InputCSR* input);


/**
 * @brief Read a CSR data from a file
 * 
 * @param filename 
 * @param verbose 
 * @return InputCSR* 
 */
InputCSR* input_read(char* filename, int verbose);

#endif