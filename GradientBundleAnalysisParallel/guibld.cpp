/***************************************************
 *                                                 *
 *  NOTE!  This file is automatically generated by *
 *         the Tecplot GUI Builder.  It is highly  *
 *         recommended that you never edit this    *
 *         file directly!                          *
 *                                                 *
 ***************************************************/



/**
 */
void BuildTab1_1(LgIndex_t  Parent)
{
  if (Tab1_1Manager != BADDIALOGID)
    return;

  Tab1_1Manager = TecGUITabAddPage(Parent,
                               "Setup");
  TecGUIFrameAdd(Tab1_1Manager,
             185,
             1663,
             9965,
             935,
                    " Sphere Mesh Parameters ");

  TecGUIFrameAdd(Tab1_1Manager,
             185,
             1299,
             10006,
             268,
                    " System Boundary Parameters ");

  TecGUIFrameAdd(Tab1_1Manager,
             2599,
             77,
             7592,
             1135,
                    " Integration Variables ");

  TecGUIFrameAdd(Tab1_1Manager,
             144,
             77,
             2331,
             1135,
                    " Select CP(s) ");

  BTNRun_BTN_T1_1 = TecGUIButtonAdd(Tab1_1Manager,
                                  8748,
                                  2651,
                                  1031,
                                  233,
                    "       Run\n",
                    BTNRun_BTN_T1_1_CB);

  TFSTPts_TF_T1_1 = TecGUITextFieldAdd(Tab1_1Manager,
                                     9346,
                                     1369,
                                     763,
                                     147,
                       TFSTPts_TF_T1_1_CB);

  LBLSTPts_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                   7138,
                                   1412,
                    "# Gradient path points:");

  LBLNumTri_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                    7035,
                                    1767,
                    "--> Num Triangles");

  TFSdRad_TF_T1_1 = TecGUITextFieldAdd(Tab1_1Manager,
                                     2104,
                                     1915,
                                     722,
                                     147,
                       TFSdRad_TF_T1_1_CB);

  LBLRad_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                 309,
                                 1776,
                    "Radial sphere approx:");

  RBRadMode_RADIO_T1_1 = TecGUIRadioBoxAdd(Tab1_1Manager,
                                         330,
                                         2105,
                                         2867,
                                         294,
                                 "Absolute",
                                 "Fraction of min CP dist.",
                                 (char *)NULL,
                                 (char *)NULL,
                                 (char *)NULL,
                                 RBRadMode_RADIO_T1_1_CB);

  TFCutoff_TF_T1_1 = TecGUITextFieldAdd(Tab1_1Manager,
                                      5735,
                                      1360,
                                      1299,
                                      147,
                       TFCutoff_TF_T1_1_CB);

  LBLCutoff_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                    3651,
                                    1412,
                    "System Rho Cutoff:");

  MLSelVars_MLST_T1_1 = TecGUIListAdd(Tab1_1Manager,
                                    2661,
                                    329,
                                    7427,
                                    875,
                       1,
                       MLSelVars_MLST_T1_1_CB);

  SCPrecise_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                   5178,
                                   173,
                                   1196,
                                   103,
                    0,
                    100,
                    0,
                    SCPrecise_SC_T1_1_CB,
                    SCPrecise_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCPrecise_SC_T1_1,FALSE);

  LBLPreci1_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                    4002,
                                    164,
                    "Precision:");

  MLSelCPs_MLST_T1_1 = TecGUIListAdd(Tab1_1Manager,
                                   226,
                                   138,
                                   2166,
                                   1048,
                       1,
                       MLSelCPs_MLST_T1_1_CB);

  LBLPrecise_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                     6953,
                                     173,
                    "Label");

  TGLsGP_TOG_T1_1 = TecGUIToggleAdd(Tab1_1Manager,
                                  412,
                                  1386,
                                  1650,
                                  112,
                               "Save GPs",
                               TGLsGP_TOG_T1_1_CB);

  TGLsGB_TOG_T1_1 = TecGUIToggleAdd(Tab1_1Manager,
                                  2083,
                                  1386,
                                  1671,
                                  112,
                               "Save GBs",
                               TGLsGB_TOG_T1_1_CB);

  LBLLevel_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                   3115,
                                   1767,
                    "Minimum number of gradient bundles (GB)");

  LBLBPGBa_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                   3115,
                                   2270,
                    "Number of bond path coincident GBs:");

  LBLGBMaxSD_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                     3115,
                                     2435,
                    "Maximum GB surface gradient path spacing:");

  LBLSDtightle_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                       3115,
                                       2105,
                    "Subdivision tightness around bond/ring features");

  TGLNoSphInt_TOG_T1_1 = TecGUIToggleAdd(Tab1_1Manager,
                                       7943,
                                       164,
                                       2909,
                                       112,
                               "Radial sphere approx.",
                               TGLNoSphInt_TOG_T1_1_CB);

  SCSDtight_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                   8830,
                                   2088,
                                   1196,
                                   103,
                    0,
                    100,
                    0,
                    SCSDtight_SC_T1_1_CB,
                    SCSDtight_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCSDtight_SC_T1_1,FALSE);

  SCBPGBs_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                 8830,
                                 2261,
                                 1196,
                                 103,
                    0,
                    100,
                    0,
                    SCBPGBs_SC_T1_1_CB,
                    SCBPGBs_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCBPGBs_SC_T1_1,FALSE);

  SCEGPDist_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                   8830,
                                   2435,
                                   1196,
                                   103,
                    0,
                    100,
                    0,
                    SCEGPDist_SC_T1_1_CB,
                    SCEGPDist_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCEGPDist_SC_T1_1,FALSE);

  LBLSDtight_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                     8170,
                                     2097,
                    "Label");

  LBLBPGBs_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                   8170,
                                   2261,
                    "Label");

  LBLEGPDistTy_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                       8170,
                                       2435,
                    "Label");

  LBLBPGBi_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                   3115,
                                   1941,
                    "Maximum GB subdivision level (around bond/ring features):");

  LBLBPGBInit_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                      8170,
                                      1923,
                    "Label");

  SCMinGBs_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                  8830,
                                  1750,
                                  1196,
                                  103,
                    0,
                    100,
                    0,
                    SCMinGBs_SC_T1_1_CB,
                    SCMinGBs_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCMinGBs_SC_T1_1,FALSE);

  SCBPGBInit_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                    8830,
                                    1915,
                                    1196,
                                    103,
                    0,
                    100,
                    0,
                    SCBPGBInit_SC_T1_1_CB,
                    SCBPGBInit_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCBPGBInit_SC_T1_1,FALSE);

  TGLSphTest_TOG_T1_1 = TecGUIToggleAdd(Tab1_1Manager,
                                      330,
                                      2452,
                                      577,
                                      112,
                               ".",
                               TGLSphTest_TOG_T1_1_CB);

  TGL_TOG_T1_1 = TecGUIToggleAdd(Tab1_1Manager,
                               495,
                               2755,
                               3961,
                               112,
                               "Enter atomic reference energies",
                               TGL_TOG_T1_1_CB);

  LBL34_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                309,
                                1949,
                    "GP seed radius:");

  TFRad_TF_T1_1 = TecGUITextFieldAdd(Tab1_1Manager,
                                   2104,
                                   1759,
                                   722,
                                   147,
                       TFRad_TF_T1_1_CB);

}

/**
 */
void BuildTab3_1(LgIndex_t  Parent)
{
  if (Tab3_1Manager != BADDIALOGID)
    return;

  Tab3_1Manager = TecGUITabAddPage(Parent,
                               "Results");
  TecGUIFrameAdd(Tab3_1Manager,
             330,
             216,
             3156,
             1213,
                    " CP Sphere ");

  TecGUIFrameAdd(Tab3_1Manager,
             392,
             1490,
             3094,
             1169,
                    "Properties");

  TecGUIFrameAdd(Tab3_1Manager,
             3631,
             207,
             6375,
             1776,
                    "Gradient Bundles");

  TecGUIFrameAdd(Tab3_1Manager,
             7118,
             2071,
             2888,
             580,
                    " Export ");

  TecGUIFrameAdd(Tab3_1Manager,
             3631,
             2071,
             3342,
             580,
                    " GB Region ");

  SLSelSphere_SLST_T3_1 = TecGUIListAdd(Tab3_1Manager,
                                      1774,
                                      251,
                                      1650,
                                      779,
                       0,
                       SLSelSphere_SLST_T3_1_CB);

  TGLSphereVis_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                        392,
                                        363,
                                        1341,
                                        112,
                               "Visible",
                               TGLSphereVis_TOG_T3_1_CB);

  BTNSphereDel_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                        412,
                                        736,
                                        845,
                                        147,
                    "Delete",
                    BTNSphereDel_BTN_T3_1_CB);

  SLSelVar_SLST_T3_1 = TecGUIListAdd(Tab3_1Manager,
                                   453,
                                   1585,
                                   2867,
                                   875,
                       0,
                       SLSelVar_SLST_T3_1_CB);

  MLSelGB_MLST_T3_1 = TecGUIListAdd(Tab3_1Manager,
                                  3693,
                                  303,
                                  3445,
                                  1594,
                       1,
                       MLSelGB_MLST_T3_1_CB);

  BTNAllGB_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                    7324,
                                    797,
                                    1877,
                                    147,
                    "         Activate All",
                    BTNAllGB_BTN_T3_1_CB);

  BTNTogMode_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                      7324,
                                      597,
                                      1877,
                                      147,
                    "      Toggle Mode",
                    BTNTogMode_BTN_T3_1_CB);

  TFNumContours_TF_T3_1 = TecGUITextFieldAdd(Tab3_1Manager,
                                           330,
                                           2746,
                                           1237,
                                           147,
                       TFNumContours_TF_T3_1_CB);

  LBL9_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                               1753,
                               2772,
                    "Number of contours");

  RBLogLin_RADIO_T3_1 = TecGUIRadioBoxAdd(Tab3_1Manager,
                                        4002,
                                        2695,
                                        1031,
                                        294,
                                 "log",
                                 "linear",
                                 (char *)NULL,
                                 (char *)NULL,
                                 (char *)NULL,
                                 RBLogLin_RADIO_T3_1_CB);

  RBCntSrc_RADIO_T3_1 = TecGUIRadioBoxAdd(Tab3_1Manager,
                                        7695,
                                        2703,
                                        1320,
                                        294,
                                 "selected",
                                 "all",
                                 (char *)NULL,
                                 (char *)NULL,
                                 (char *)NULL,
                                 RBCntSrc_RADIO_T3_1_CB);

  LBLCntSrc_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                    5261,
                                    2764,
                    "Make contour values from");

  LBL13_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                8830,
                                2764,
                    "sphere(s)");

  BTNExport_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                     7427,
                                     2443,
                                     1650,
                                     147,
                    "Export to CSV",
                    BTNExport_BTN_T3_1_CB);

  TGLExGBs_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                    7303,
                                    2149,
                                    3424,
                                    112,
                               "Include ALL individual GBs",
                               TGLExGBs_TOG_T3_1_CB);

  TGLShowMesh_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                       392,
                                       537,
                                       1836,
                                       112,
                               "Show Mesh",
                               TGLShowMesh_TOG_T3_1_CB);

  BTNSelGB_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                    5302,
                                    2426,
                                    1217,
                                    147,
                    "Select GB",
                    BTNSelGB_BTN_T3_1_CB);

  TFGrpNum_TF_T3_1 = TecGUITextFieldAdd(Tab3_1Manager,
                                      4477,
                                      2426,
                                      722,
                                      147,
                       TFGrpNum_TF_T3_1_CB);

  LBL20_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                3693,
                                2469,
                    "Group #");

  LBL21_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                3693,
                                2131,
                    "1. Activate boundary w/ toggle mode\n2. Select GB within boundary");

  BTNFndBas_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                     7303,
                                     311,
                                     2001,
                                     233,
                    "Find Special \nGradient Bundles",
                    BTNFndBas_BTN_T3_1_CB);

  BTNSmooth_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                     2331,
                                     2487,
                                     990,
                                     147,
                    "Smooth",
                    BTNSmooth_BTN_T3_1_CB);

  TGLSmoothAll_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                        495,
                                        2513,
                                        1815,
                                        112,
                               "All spheres",
                               TGLSmoothAll_TOG_T3_1_CB);

  SCRad_SC_T3_1 = TecGUIScaleAdd(Tab3_1Manager,
                               453,
                               1256,
                               2599,
                               103,
                    0,
                    100,
                    0,
                    SCRad_SC_T3_1_CB,
                    SCRad_SCD_T3_1_CB);


  TecGUIScaleShowNumericDisplay(SCRad_SC_T3_1,FALSE);

  TGLRadAll_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                     474,
                                     1091,
                                     1836,
                                     112,
                               "All Spheres",
                               TGLRadAll_TOG_T3_1_CB);

  TGLRadAbs_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                     1898,
                                     1091,
                                     1568,
                                     112,
                               "Absolute",
                               TGLRadAbs_TOG_T3_1_CB);

  LBL27_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                474,
                                918,
                    "Radius:");

  LBLRadLab_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                    3094,
                                    1247,
                    "Label");

  LBL29_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                7324,
                                2305,
                    "Uncheck to only include active.");

  TGL_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                               7365,
                               1057,
                               536,
                               112,
                               ".",
                               TGL_TOG_T3_1_CB);

  LBL31_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                7654,
                                996,
                    "Save special gradient\nbundle surfaces");

}

/**
 */
void BuildDialog1(LgIndex_t  ParentDialog)
{
  if (Dialog1Manager != BADDIALOGID)
    return;

  Dialog1Manager = TecGUIDialogCreateModeless(ParentDialog,
                                             10550,
                                             3262,
                                             "Gradient Bundle Analysis",
                                             Dialog1Init_CB,
                                             Dialog1CloseButton_CB,
                                             NULL)
;  TAB1_TB_D1 = TecGUITabAdd(Dialog1Manager,
                          41,
                          181,
                          10460,
                          3050,
                    TAB1_TBA_D1_CB,
                    TAB1_TBD_D1_CB);

  BuildTab1_1(TAB1_TB_D1);
  BuildTab3_1(TAB1_TB_D1);
}


/**
 */
void InitTGB(void)
{
/* Currently not used */
}
