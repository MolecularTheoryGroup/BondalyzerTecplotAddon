#ifndef STREAMTRACEOPERATIONS_H_
#define STREAMTRACEOPERATIONS_H_



Boolean_t StreamtraceResample(EntIndex_t STZoneNum, LgIndex_t NumSTPoints);
Boolean_t StreamtraceResampleToExistingZone(EntIndex_t STZoneNum,
	EntIndex_t NewZoneNum,
	LgIndex_t NumSTPoints,
	EntIndex_t NumVars,
	EntIndex_t* XYZVarNums,
	vec3 & MinVolPt);

Boolean_t StreamtraceConcatenateResample(EntIndex_t STZoneNum, EntIndex_t VolGPNum, LgIndex_t NumSTPoints);
Boolean_t StreamtraceConcatenateResampleToExistingZone(EntIndex_t STZoneNum,
	EntIndex_t NewZoneNum,
	EntIndex_t VolGPNum,
	LgIndex_t NumSTPoints,
	EntIndex_t NumVars,
	EntIndex_t* XYZVarNums,
	vec3 & MinVolPt);

#endif