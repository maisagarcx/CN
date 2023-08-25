#include <string.h>  
 #include <stdio.h>  
 #include <stdlib.h>  
 #define MAX 100  
 void readSizeOfMatrix(int *R, int *C){
     do{
         printf("Please, input the how much rows do you want: ");          
     scanf("%d", R);
     if(*R>MAX){
          printf("Ops, this matrix has a MAX of 100 rows. Try again.\n");
          }  
     }while(*R>MAX);    
     do{
     printf("Please, input the how much columns do you want: ");          
     scanf("%d", C);
     if(*C>MAX){
          printf("Ops, this matrix has a MAX of 100 columns. Try again.\n");
          }  
     }while(*R>MAX);      
 }  
 void readElementsOfMatrix(int R, int C, int M[MAX][MAX]){     
     for(int i=0;i<R;i++){ 
         for(int j=0;j<C;j++){                  
             printf("Input the element matrix[%d][%d]: ", i, j); 
                          scanf("%d", &M[i][j]); 
         }          
     }  
 }  
 void displayElementsOfMatrix(int R, int C, int M[MAX][MAX]){
     printf("Your matrix is:\n"); 
     for(int i=0;i<R;i++){ 
         for(int j=0;j<C;j++){ 
             printf("%d ", M[i][j]); 
         }                  
         printf("\n"); 
     }  
 }  
 void findHighestNumberInEveryRow(int R, int C, int M[MAX][MAX], int H[MAX]){ 
     int aux; 
     for(int i=0;i<R;i++){ 
         aux = M[i][0]; 
                  for(int j=0;j<C;j++){  
             if(M[i][j]>=aux){ 
                 aux=M[i][j]; 
                          } 
         } 
         H[i]=aux;  
     }  
 }  
 void displayArrayOfHighests(int R,int H[MAX]){  
     for(int i=0;i<R;i++){ 
         printf("The highest in the row %d is %d.\n", i, H[i]); 
     }  
 }  
 int main(){ 
     //this program will find the highest element in every row in a matrix          
     int rows, columns; 
     int matrix[MAX][MAX], highests[MAX]; 
     readSizeOfMatrix(&rows, &columns); 
     readElementsOfMatrix(rows, columns, matrix); 
     displayElementsOfMatrix(rows, columns, matrix); 
     findHighestNumberInEveryRow(rows, columns, matrix, highests); 
     displayArrayOfHighests(rows, highests);  
 }