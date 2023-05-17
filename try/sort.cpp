#include <iostream>
#include <algorithm>
#include "analyzer.hpp"



int main(int argc, char** argv){
	int arr[] = {5, 6, 4, 1, 7, 9, 2, 3, 4, 1, 8, 9};
	int size = sizeof(arr)/sizeof(int);
	int i;

	std::sort(arr, arr+size);
	
	for (i = 0; i < size; i++){
		std::cout << arr[i] << std::endl;
	}
	std::cout << std::endl;

	return 0;
}
