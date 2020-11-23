#pragma  once

typedef enum
{
    StereoVisOp_SinglePt,
    StereoVisOp_AllCps,
    StereoVisOp_PlotSingleZoneTest,
    END_StereoVisOp_e,
    StereoVisOp_Invalid = BadEnumValue
} StereoVisOp_e;

void stereoVis(StereoVisOp_e stereoVisOp);
void plotSingleZone(); // for testing
