#ifndef GTAENGINE_H_
#define GTAENGINE_H_





Boolean_t GTARunGTA(std::vector<vec3> *Pts);
Boolean_t PopulateRadiusVar(EntIndex_t RadVarNum,
	EntIndex_t XYZVarNums[3],
	EntIndex_t CutoffVarNum,
	double     CutoffValue,
	EntIndex_t FEZoneNum,
	vec3 Origin,
	vec3 Axis);
Boolean_t SetPercent(unsigned int CurrentNum, unsigned int TotalNum, const char* ProgresssText);

#endif