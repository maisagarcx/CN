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
void whichNorm(int *normPointer){
	printf("Please choose a P-norm for your calculus: ");
	do{
		scanf("%d", normPointer);
	}while(*normPointer>=LIMIT);
}
float calculateResult(int size, int norm, int vector[MAX]){
	float sum = 0.0, result;
	float exponent = 1/norm;
	for(int index=0;index<size;index++){
		sum+=pow(abs(vector[index]), norm);
	}
	return pow(sum, exponent);
}
void displayResult(int size, int norm, int array, float result){
	
}
void main(){
	//this program will calculate a norm for any
	int size, norm, array[MAX];
	readSizeOfArray(&size);
	readElementsOfArray(size, array);
	displayElementsOfArray(size, array);
	whichNorm(&norm);
	float result = calculateResult(size, norm, array);
	displayResult(size, norm, array, result);
}
