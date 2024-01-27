#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 100
void readSizeOfArray(int *size){
    do{
        printf("Please, input the size of your array: ");
        scanf("%d", size);
        if(*size>MAX){
            printf("Ops, this array has a MAX of 100 elements. Try again.\n");
        }
    }while(*size>MAX);
}
void readElementsOfArray(int size, int array[MAX]){
    printf("Please, input the elements of your array: ");
    for(int i=0;i<size;i++){
        printf("Array[%d]: ", i);
        scanf("%d", &array[i]);
    }
}
void displayElementsOfArray(int size, int array[MAX]){
    printf("Your array is: ");
    for(int i=0;i<size;i++){
        printf("%d ", array[i]);
    }
}
void calculateAndDisplayTheAverage(int size, int array[MAX], float *average){
    float sum = 0;
    for(int i=0;i<size;i++){
        sum+=array[i];
    }
    *average=sum/size;
    printf("Your average is %f", *average);
}
void calculateAndDisplayTheStandardDeviation(int size, int array[MAX], float average, float *standardDeviation){
    float sum = 0;
    for(int i=0;i<size;i++){
        sum+=pow(array[i], 2);
    }
    *standardDeviation=sqrt((sum-(pow(average, 2)/size))/(size - 1));
    printf("Your standard deviation is %f", *standardDeviation);
}
void main() {
    //this code calculates the average and the standard deviation of a array
    int size, array[MAX];
    float average, standardDeviation;
    readSizeOfArray(&size);
    readElementsOfArray(size, array);
    displayElementsOfArray(size, array);
    calculateAndDisplayTheAverage(size, array, &average);
    calculateAndDisplayTheStandardDeviation(size, array, average, &standardDeviation);
}
