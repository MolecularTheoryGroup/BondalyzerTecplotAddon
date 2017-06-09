/*
 * main.cpp
 *
 *  Created on: Aug 8, 2015
 *      Author: Haiiro
 */

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(){

	ifstream Input("/Users/Haiiro/Desktop/cubane.txt");

	vector<double> Vals;
	Vals.reserve(20000);

	double Tmp;

	while (!Input.eof()){
		Input >> Tmp;
		Vals.push_back(Tmp);
	}

	Input.close();

	ofstream Output("/Users/Haiiro/Desktop/out.csv", ios::out | ios::trunc);

	int NumCols = 20;
	int ColNum = 0;

	for (double & i : Vals){
		ColNum++;
		if (ColNum >= NumCols){
			Output << i << endl;
			ColNum = 0;
		}
		else
			Output << i << ",";
	}

	Output.close();

	return 0;
}

