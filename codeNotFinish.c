#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 1024
#define LIMIT 0
void readSizeOfArray(int *sizePointer){ 
    printf("Please, input the size of your array: ");
    do{
    	scanf("%d", sizePointer); 
	}while(*sizePointer<=LIMIT||*sizePointer>MAX);
}
void readElementsOfArray(int yourSize, int vector[MAX]){ 
    for(int index=0;index<yourSize;index++){ 
        printf("Input the element vector[%d]: ", index); 
        scanf("%d", &vector[index]); 
    }  
} 
void displayElementsOfArray(int yourSize, int vector[MAX]){ 
    printf("\nYour array is:\n"); 
    for(int index=0;index<yourSize;index++){ 
        printf("%d ", vector[index]); 
    } 
}
void whichNorm(*normPointer){
	printf("Please choose a P-norm for your calculus: ");
	do{
		scanf("%d", &normPointer);
	}while(*normPointer>=LIMIT);
}
void 
void main(){
	//this program will calculate a norm for any
	int size, norm, array[MAX];
	readSizeOfArray(&size);
	readElementsOfArray(size, array);
	displayElementsOfArray(size, array);
	whichNorm(&norm);
	sumOfAny(size, norm, array);
	
}
