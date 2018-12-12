#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include <string>
#include <vector>

using std::vector;
using std::string;

static int IntPrecise = 2;
static vector<string> IntPrecisionLabels = { "Coarse", "Normal", "Fine", "Finer", "Finest" };

// void GBAIntegrationPrepareGUI();

int GetRecommendedIntPrecision();

Boolean_t PerformIntegration(vector<string> const & AtomNameList,
	vector<string> const & IntVarNameList,
	vector<int> const & IntVarNumList,
	Boolean_t IntegrateVolume,
	int IntResolution,
	Boolean_t ActiveGBsOnly);

#endif