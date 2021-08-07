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
             196,
             1662,
             9970,
             939,
                    " Sphere Mesh Parameters ");

  TecGUIFrameAdd(Tab1_1Manager,
             196,
             1300,
             10000,
             266,
                    " System Boundary Parameters ");

  TecGUIFrameAdd(Tab1_1Manager,
             2583,
             76,
             7598,
             1135,
                    " Integration Variables ");

  TecGUIFrameAdd(Tab1_1Manager,
             135,
             82,
             2326,
             1135,
                    " Select CP(s) ");

  BTNRun_BTN_T1_1 = TecGUIButtonAdd(Tab1_1Manager,
                                  8731,
                                  2652,
                                  1027,
                                  234,
                    "       Run\n",
                    BTNRun_BTN_T1_1_CB);

  TFSTPts_TF_T1_1 = TecGUITextFieldAdd(Tab1_1Manager,
                                     9350,
                                     1370,
                                     770,
                                     145,
                       TFSTPts_TF_T1_1_CB);

  LBLSTPts_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                   7145,
                                   1414,
                    "# Gradient path points:");

  LBLNumTri_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                    7039,
                                    1770,
                    "--> Num Triangles");

  TFRad_TF_T1_1 = TecGUITextFieldAdd(Tab1_1Manager,
                                   1268,
                                   1833,
                                   710,
                                   145,
                       TFRad_TF_T1_1_CB);

  LBLRad_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                 377,
                                 1852,
                    "Radius:");

  RBRadMode_RADIO_T1_1 = TecGUIRadioBoxAdd(Tab1_1Manager,
                                         407,
                                         2055,
                                         2870,
                                         298,
                                 "Absolute",
                                 "Fraction of min CP dist.",
                                 (char *)NULL,
                                 (char *)NULL,
                                 (char *)NULL,
                                 RBRadMode_RADIO_T1_1_CB);

  TFCutoff_TF_T1_1 = TecGUITextFieldAdd(Tab1_1Manager,
                                      5725,
                                      1364,
                                      1299,
                                      145,
                       TFCutoff_TF_T1_1_CB);

  LBLCutoff_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                    3655,
                                    1408,
                    "System Rho Cutoff:");

  MLSelVars_MLST_T1_1 = TecGUIListAdd(Tab1_1Manager,
                                    2673,
                                    323,
                                    7432,
                                    869,
                       1,
                       MLSelVars_MLST_T1_1_CB);

  SCPrecise_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                   5181,
                                   177,
                                   1193,
                                   107,
                    0,
                    100,
                    0,
                    SCPrecise_SC_T1_1_CB,
                    SCPrecise_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCPrecise_SC_T1_1,FALSE);

  LBLPreci1_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                    4003,
                                    171,
                    "Precision:");

  MLSelCPs_MLST_T1_1 = TecGUIListAdd(Tab1_1Manager,
                                   211,
                                   139,
                                   2160,
                                   1053,
                       1,
                       MLSelCPs_MLST_T1_1_CB);

  LBLPrecise_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                     6949,
                                     177,
                    "Label");

  TGLsGP_TOG_T1_1 = TecGUIToggleAdd(Tab1_1Manager,
                                  407,
                                  1383,
                                  1722,
                                  126,
                               "Save GPs",
                               TGLsGP_TOG_T1_1_CB);

  TGLsGB_TOG_T1_1 = TecGUIToggleAdd(Tab1_1Manager,
                                  2084,
                                  1383,
                                  1737,
                                  126,
                               "Save GBs",
                               TGLsGB_TOG_T1_1_CB);

  LBLLevel_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                   2734,
                                   1770,
                    "Minimum number of gradient bundles (GB)");

  LBLBPGBa_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                   2734,
                                   2271,
                    "Number of bond path coincident GBs:");

  LBLGBMaxSD_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                     2734,
                                     2436,
                    "Maximum GB surface gradient path spacing:");

  LBLSDtightle_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                       2734,
                                       2100,
                    "Subdivision tightness around bond/ring features");

  TGLNoSphInt_TOG_T1_1 = TecGUIToggleAdd(Tab1_1Manager,
                                       7946,
                                       158,
                                       2991,
                                       126,
                               "Radial sphere approx.",
                               TGLNoSphInt_TOG_T1_1_CB);

  SCSDtight_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                   8822,
                                   2093,
                                   1193,
                                   107,
                    0,
                    100,
                    0,
                    SCSDtight_SC_T1_1_CB,
                    SCSDtight_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCSDtight_SC_T1_1,FALSE);

  SCBPGBs_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                 8822,
                                 2265,
                                 1193,
                                 107,
                    0,
                    100,
                    0,
                    SCBPGBs_SC_T1_1_CB,
                    SCBPGBs_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCBPGBs_SC_T1_1,FALSE);

  SCEGPDist_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                   8822,
                                   2436,
                                   1193,
                                   107,
                    0,
                    100,
                    0,
                    SCEGPDist_SC_T1_1_CB,
                    SCEGPDist_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCEGPDist_SC_T1_1,FALSE);

  LBLSDtight_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                     8172,
                                     2093,
                    "Label");

  LBLBPGBs_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                   8172,
                                   2265,
                    "Label");

  LBLEGPDistTy_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                       8172,
                                       2436,
                    "Label");

  LBLBPGBi_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                   2734,
                                   1935,
                    "Maximum GB subdivision level (around bond/ring features):");

  LBLBPGBInit_LBL_T1_1 = TecGUILabelAdd(Tab1_1Manager,
                                      8172,
                                      1928,
                    "Label");

  SCMinGBs_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                  8822,
                                  1744,
                                  1193,
                                  107,
                    0,
                    100,
                    0,
                    SCMinGBs_SC_T1_1_CB,
                    SCMinGBs_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCMinGBs_SC_T1_1,FALSE);

  SCBPGBInit_SC_T1_1 = TecGUIScaleAdd(Tab1_1Manager,
                                    8822,
                                    1916,
                                    1193,
                                    107,
                    0,
                    100,
                    0,
                    SCBPGBInit_SC_T1_1_CB,
                    SCBPGBInit_SCD_T1_1_CB);


  TecGUIScaleShowNumericDisplay(SCBPGBInit_SC_T1_1,FALSE);

  TGLSphTest_TOG_T1_1 = TecGUIToggleAdd(Tab1_1Manager,
                                      392,
                                      2404,
                                      664,
                                      126,
                               ".",
                               TGLSphTest_TOG_T1_1_CB);

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
             332,
             215,
             3142,
             1211,
                    " CP Sphere ");

  TecGUIFrameAdd(Tab3_1Manager,
             392,
             1491,
             3096,
             1167,
                    "Properties");

  TecGUIFrameAdd(Tab3_1Manager,
             3640,
             203,
             6374,
             1782,
                    "Gradient Bundles");

  TecGUIFrameAdd(Tab3_1Manager,
             7130,
             2068,
             2885,
             583,
                    " Export ");

  TecGUIFrameAdd(Tab3_1Manager,
             3640,
             2074,
             3353,
             577,
                    " GB Region ");

  SLSelSphere_SLST_T3_1 = TecGUIListAdd(Tab3_1Manager,
                                      1767,
                                      253,
                                      1646,
                                      780,
                       0,
                       SLSelSphere_SLST_T3_1_CB);

  TGLSphereVis_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                        392,
                                        361,
                                        1420,
                                        126,
                               "Visible",
                               TGLSphereVis_TOG_T3_1_CB);

  BTNSphereDel_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                        422,
                                        742,
                                        845,
                                        145,
                    "Delete",
                    BTNSphereDel_BTN_T3_1_CB);

  SLSelVar_SLST_T3_1 = TecGUIListAdd(Tab3_1Manager,
                                   468,
                                   1586,
                                   2885,
                                   869,
                       0,
                       SLSelVar_SLST_T3_1_CB);

  MLSelGB_MLST_T3_1 = TecGUIListAdd(Tab3_1Manager,
                                  3701,
                                  304,
                                  3444,
                                  1598,
                       1,
                       MLSelGB_MLST_T3_1_CB);

  BTNAllGB_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                    7311,
                                    793,
                                    1888,
                                    145,
                    "         Activate All",
                    BTNAllGB_BTN_T3_1_CB);

  BTNTogMode_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                      7311,
                                      596,
                                      1888,
                                      145,
                    "      Toggle Mode",
                    BTNTogMode_BTN_T3_1_CB);

  TFNumContours_TF_T3_1 = TecGUITextFieldAdd(Tab3_1Manager,
                                           332,
                                           2747,
                                           1238,
                                           145,
                       TFNumContours_TF_T3_1_CB);

  LBL9_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                               1752,
                               2772,
                    "Number of contours");

  RBLogLin_RADIO_T3_1 = TecGUIRadioBoxAdd(Tab3_1Manager,
                                        4003,
                                        2696,
                                        1027,
                                        298,
                                 "log",
                                 "linear",
                                 (char *)NULL,
                                 (char *)NULL,
                                 (char *)NULL,
                                 RBLogLin_RADIO_T3_1_CB);

  RBCntSrc_RADIO_T3_1 = TecGUIRadioBoxAdd(Tab3_1Manager,
                                        7689,
                                        2702,
                                        1329,
                                        298,
                                 "selected",
                                 "all",
                                 (char *)NULL,
                                 (char *)NULL,
                                 (char *)NULL,
                                 RBCntSrc_RADIO_T3_1_CB);

  LBLCntSrc_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                    5257,
                                    2766,
                    "Make contour values from");

  LBL13_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                8837,
                                2766,
                    "sphere(s)");

  BTNExport_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                     7417,
                                     2449,
                                     1646,
                                     145,
                    "Export to CSV",
                    BTNExport_BTN_T3_1_CB);

  TGLExGBs_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                    7296,
                                    2144,
                                    3519,
                                    126,
                               "Include ALL individual GBs",
                               TGLExGBs_TOG_T3_1_CB);

  TGLShowMesh_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                       392,
                                       532,
                                       1933,
                                       126,
                               "Show Mesh",
                               TGLShowMesh_TOG_T3_1_CB);

  BTNSelGB_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                    5302,
                                    2430,
                                    1223,
                                    145,
                    "Select GB",
                    BTNSelGB_BTN_T3_1_CB);

  TFGrpNum_TF_T3_1 = TecGUITextFieldAdd(Tab3_1Manager,
                                      4471,
                                      2430,
                                      710,
                                      145,
                       TFGrpNum_TF_T3_1_CB);

  LBL20_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                3701,
                                2468,
                    "Group #");

  LBL21_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                3701,
                                2131,
                    "1. Activate boundary w/ toggle mode\n2. Select GB within boundary");

  BTNFndBas_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                     7296,
                                     310,
                                     1994,
                                     234,
                    "Find Special \nGradient Bundles",
                    BTNFndBas_BTN_T3_1_CB);

  BTNSmooth_BTN_T3_1 = TecGUIButtonAdd(Tab3_1Manager,
                                     2326,
                                     2487,
                                     997,
                                     145,
                    "Smooth",
                    BTNSmooth_BTN_T3_1_CB);

  TGLSmoothAll_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                        498,
                                        2506,
                                        1888,
                                        126,
                               "All spheres",
                               TGLSmoothAll_TOG_T3_1_CB);

  SCRad_SC_T3_1 = TecGUIScaleAdd(Tab3_1Manager,
                               453,
                               1256,
                               2598,
                               107,
                    0,
                    100,
                    0,
                    SCRad_SC_T3_1_CB,
                    SCRad_SCD_T3_1_CB);


  TecGUIScaleShowNumericDisplay(SCRad_SC_T3_1,FALSE);

  TGLRadAll_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                     483,
                                     1084,
                                     1903,
                                     126,
                               "All Spheres",
                               TGLRadAll_TOG_T3_1_CB);

  TGLRadAbs_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                                     1888,
                                     1084,
                                     1646,
                                     126,
                               "Absolute",
                               TGLRadAbs_TOG_T3_1_CB);

  LBL27_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                483,
                                919,
                    "Radius:");

  LBLRadLab_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                    3096,
                                    1249,
                    "Label");

  LBL29_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                7326,
                                2303,
                    "Uncheck to only include active.");

  TGL_TOG_T3_1 = TecGUIToggleAdd(Tab3_1Manager,
                               7372,
                               1046,
                               604,
                               126,
                               ".",
                               TGL_TOG_T3_1_CB);

  LBL31_LBL_T3_1 = TecGUILabelAdd(Tab3_1Manager,
                                7643,
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
                          45,
                          183,
                          10453,
                          3051,
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
