
#include <string>
#include <vector>

#include "TECADDON.h"
#include "ZONEVARINFO.h"

using std::vector;
using std::string;

void LogDataValues(vector<double> &LevelVector, double Min, double Max, double BaseValue)
{
	if (Max - Min <= 0 && Min != 1e10)
	{
		LevelVector.push_back(0.0);
	}
	else
	{
		int LevelInt = 0, MultInt = 1, ExpInt = 0;
		double TempValue = 0;
		Boolean_t BreakLoop = FALSE;
		TempValue = 1.0;
		if (Min != 1e10)
		{
			while (!BreakLoop)
			{
				if (TempValue > Min) TempValue *= 0.1;
				else if (TempValue * 10 <= Min) TempValue *= 10;
				else
				{
					BaseValue = TempValue;
					BreakLoop = TRUE;
				}
			}
		}
		BreakLoop = FALSE;
		while (!BreakLoop)
		{
			TempValue = MultInt * BaseValue * pow(10.0, ExpInt);
			if (TempValue >= Max)
			{
				TempValue = Max;
				BreakLoop = TRUE;
			}
			LevelVector.push_back(TempValue);

			if (MultInt >= 9)
			{
				MultInt = 1;
				ExpInt++;
			}
			else MultInt++;
		}
	}
}

void X2DataValues(vector<double> &LevelVector, double Min, double Max)
{
	if (Max - Min <= 0 && Min != 1e10)
	{
		LevelVector.push_back(0.0);
	}
	else
	{
		LevelVector.push_back(Min);
		while (LevelVector.back() < Max){
			LevelVector.push_back(2 * LevelVector.back());
		}
	}
}

string ReplaceStringChar(string const & InString, 
	string const & OldChar, 
	string const & NewChar)
{
	string NewString = InString;

	size_t Pos = NewString.find_first_of(OldChar);
	while (Pos != string::npos){
		NewString.erase(Pos, OldChar.length());
		NewString.insert(Pos, NewChar);
		Pos = NewString.find_first_of(OldChar);
	}

	return NewString;
}

string RemoveStringChar(string const & InString,
	string const & OldChar)
{
	string NewString = InString;

	size_t Pos = NewString.find_first_of(OldChar);
	while (Pos != string::npos){
		NewString.erase(Pos, OldChar.length());
		Pos = NewString.find_first_of(OldChar);
	}

	return NewString;
}