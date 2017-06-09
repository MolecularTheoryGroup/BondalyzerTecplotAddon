#ifndef GUICB_H_
#define GUICB_H_



void GBAResultViewerPrepareGUI();
Boolean_t GBAProcessSystemPrepareGUI();
void GBAProcessSystemUpdateNumTriangles();
void GBAProcessSystemUpdateNumGPsPerGB();
void GBAProcessSystemLabelSelectedCPs();
void GBAProcessSystemDeleteCPLabels();

void PrepareIntegration(Boolean_t IntegratingFromIntTab);

void GBAReloadDialog();

#endif