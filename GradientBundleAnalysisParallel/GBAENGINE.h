#pragma once

#include <string>
#include <vector>
#include <map>

#include "CSM_FE_VOLUME.h"

using std::vector;
using std::string;



void NewMainFunction();
void MainFunction();
void GradPathTest();
void CreateCircularGBsGetUserInfo();

void FindSphereBasins();

void ComputeGradientOnSphereSurface(int SphereZoneNum, vector<int> XYZVarNums, int ValVarNum);