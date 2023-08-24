#include <string.h> 
#include <stdio.h> 
#include <stdlib.h> 
#define MAX 100 
void readSize(int *L, int *C){ 	
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
void findHighestNumberInEveryLine(int L, int C, int M[MAX][MAX], int H[MAX]){
    int aux;
    for(int i=0;i<L;i++){
        aux = M[0][i];
 		for(int j=0;j<C;j++){ 
            if(M[i][j]>=aux){
                aux=M[i][j];
 			}
        }
        H[i]=aux; 
    } 
} 
void showHighests(int H[MAX]){ 
    for(int i=0;H[i]!='\0';i++){
        printf("The highest in the line %d is %d.\n", i, H[i]);
 	} 
} int main(){
    //this program will find the highest element in every line in a matrix 	
    int lines, columns;
    int matrix[MAX][MAX], highests[MAX];
    readSize(&lines, &columns);
    readMatrix(lines, columns, matrix);
    showMatrix(lines, columns, matrix);
    findHighestNumberInEveryLine(lines, columns, matrix, highests);
    showHighests(highests); 
}

