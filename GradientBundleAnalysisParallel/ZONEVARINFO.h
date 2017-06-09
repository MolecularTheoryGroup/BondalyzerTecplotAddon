#ifndef ZONEVARINFO_H_
#define ZONEVARINFO_H_

#include <vector>
#include <string>

using std::vector;
using std::string;

void LogDataValues(vector<double> &LevelVector, double Min, double Max, double BaseValue);
void X2DataValues(vector<double> &LevelVector, double Min, double Max);

const string ReplaceStringChar(const string & InString,
	const string & OldChar,
	const string & NewChar);
const string RemoveStringChar(const string & InString,
	const string & OldChar);

#endif