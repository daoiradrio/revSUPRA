#include <iostream>
#include <algorithm>



#define TOL 0.0001



void dac(int* arr, int search, int l, int r){
	int m;
	while (l <= r){
		m = l + (r - l)/2;
		if (search == arr[m]){
			std::cout << "Found!" << std::endl;
			return;
		}
		else if (search < arr[m]){
			r = m - 1;
		}
		else{
			l = m + 1;
		}
	}
	std::cout << "Not found..." << std::endl;
	return;
}



int main(int argc, char** argv){
	//int arr[] = {5, 6, 42, 10, 7, 9, 2, 3, 4, 1, 8, 9, 7}; // ungerade Anzahl, keine Duplikate
	//int arr[] = {5, 6, 4, 1, 7, 9, 2, 3, 4, 1, 8, 9, 7}; // ungerade Anzahl, Duplikate
	//int arr[] = {5, 6, 42, 10, 7, 9, 2, 3, 4, 1, 8, 9}; // gerade Anzahl, keine Duplikate
	int arr[] = {5, 6, 4, 1, 7, 9, 2, 3, 4, 1, 8, 9}; // gerade Anzahl, Duplikate
	int n = sizeof(arr)/sizeof(int);
	int i, j;

	std::sort(arr, arr+n);
	
	for (i = 0; i < n; i++){
		std::cout << arr[i] << std::endl;
	}
	std::cout << std::endl;

	dac(arr, 50, 0, n);

	return 0;
}
