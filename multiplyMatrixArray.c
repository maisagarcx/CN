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
 void readSizeOfArray(int *S, int L){ 
         printf("Please, input the size of your array: "); 
     do{ 
         scanf("%d", S); 
         if(*S!=L){ 
             printf("Your array must have the same number of lines in your matrix.\n"); 
             printf("PS: Your matrix have %d lines.\n", L);
             printf("Try again."); 
         } 
     }while(*S!=L); 
 }  
 void readElementsOfArray(int S, int V[MAX]){ 
     for(int i=0;i<S;i++){ 
         printf("Input the element vector[%d]: ", i); 
         scanf("%d", &V[i]); 
     }  
 } 
 void displayElementsOfArray(int S, int V[MAX]){ 
     printf("\nYour array is:\n"); 
     for(int i=0;i<S;i++){ 
         printf("%d ", V[i]); 
     } 
 } 
 void multiplyMatrixAndArray(int R, int C, int M[MAX][MAX], int V[MAX], int NV[MAX]){ 
     for(int i=0;i<R;i++){ 
         int sum = 0; 
         for(int j=0;j<C;j++){ 
             sum+=M[i][j]*V[j]; 
         } 
         NV[i]=sum; 
     } 
 } 
 int main(){ 
     //this program will multiply a matrix and a array
     int rows, columns, size; 
     int matrix[MAX][MAX], array[MAX], newArray[MAX]; 
     readSizeOfMatrix(&rows, &columns); 
     readElementsOfMatrix(rows, columns, matrix); 
     displayElementsOfMatrix(rows, columns, matrix); 
     readSizeOfArray(&size, rows); 
     readElementsOfArray(size, array); 
     displayElementsOfArray(); 
     multiplyMatrixAndArray(rows, columns, matrix, array, newArray); 
     displayElementsOfArray(size, newArray); 
 }
