#pragma once

extern bool GBA_REINIT_DIALOG;

void GBAResultViewerPrepareGUI();
Boolean_t GBAProcessSystemPrepareGUI();
void GBAProcessSystemUpdateNumTriangles();
void GBAProcessSystemLabelSelectedCPs();
void GBAProcessSystemDeleteCPLabels();

void PrepareIntegration(Boolean_t IntegratingFromIntTab);

void GBAReloadDialog(); 
