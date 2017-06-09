#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include <string>
#include <vector>

using std::vector;
using std::string;

static int IntPrecise = 1;
static vector<string> IntPrecisionLabels = { "Coarse", "Normal", "Fine", "Finer", "Finest" };

void GBAIntegrationPrepareGUI();

const int GetRecommendedIntPrecision();

const Boolean_t PerformIntegration(const vector<string> & AtomNameList,
	const vector<string> & IntVarNameList,
	const vector<int> & IntVarNumList,
	const Boolean_t & IntegrateVolume,
	const int & IntResolution);

#endif