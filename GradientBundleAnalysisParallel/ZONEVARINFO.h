#pragma once

#include <vector>
#include <string>

using std::vector;
using std::string;

void LogDataValues(vector<double> &LevelVector, double Min, double Max, double BaseValue);
void X2DataValues(vector<double> &LevelVector, double Min, double Max);

string ReplaceStringChar(string const & InString,
	string const & OldChar,
	string const & NewChar);
string RemoveStringChar(string const & InString,
	string const & OldChar);
