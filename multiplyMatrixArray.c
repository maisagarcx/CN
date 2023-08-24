#include <string.h> 
#include <stdio.h> 
#include <stdlib.h> 
#define MAX 100 

void readSizeOfMatrix(int *L, int *C){ 	
    printf("Input the how much lines do you want: "); 	
    scanf("%d", L); 	
    printf("Input the how much columns do you want: "); 	
    scanf("%d", C); 
} 
void readMatrix(int L, int C, int M[MAX][MAX]){
    getchar(); 	
    for(int i=0;i<L;i++){
        for(int j=0;j<C;j++){ 		
            printf("Input the element matrix[%d][%d]: ", i, j);
 			scanf("%d", &M[i][j]);
        } 	
    } 
} 
void showMatrix(int L, int C, int M[MAX][MAX]){
    printf("Your matrix is:\n");
    for(int i=0;i<L;i++){
        for(int j=0;j<C;j++){
            printf("%d ", M[i][j]);
        } 		
        printf("\n");
    } 
} 
void readSizeOfVector(int *S, int L){
	printf("Input the size of your vector: ");
    do{
        scanf("%d", S);
        if(*S!=L){
            printf("Your vector must have the same number of lines in your matrix.\n");
            printf("PS: Your matrix have %d lines.\n", L);
            printf("Try again.");
        }
    }while(*S!=L);
} 
void readVector(int S, int V[MAX]){
    for(int i=0;i<S;i++){
        printf("Input the element vector[%d]: ", i);
        scanf("%d", &V[i]);
    } 
}
void showVector(int V[MAX]){
    printf("\nYour vector is:\n");
    for(int i=0;V[i]!='\0';i++){
        printf("%d ", V[i]);
    }
}
void multipleMatrixAndVector(int L, int C, int M[MAX][MAX], int V[MAX], int NV[MAX]){
    for(int i=0;i<L;i++){
        int sum = 0;
        for(int j=0;j<C;j++){
            sum+=M[i][j]*V[j];
        }
        NV[i]=sum;
    }
}
int main(){
    //this program will multiple a matrix and a vector
    int lines, columns, size;
    int matrix[MAX][MAX], vector[MAX], newVector[MAX];
    readSizeOfMatrix(&lines, &columns);
    readMatrix(lines, columns, matrix);
    showMatrix(lines, columns, matrix);
    readSizeOfVector(&size, lines);
    readVector(size, vector);
    showVector(vector);
    multipleMatrixAndVector(lines, columns, matrix, vector, newVector);
    showVector(newVector);
}

