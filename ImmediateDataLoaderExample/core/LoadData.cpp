#include "tptoolbox.h"
#include "StringList.h"
#include "Set.h"
#include "../Lock.h"
#include <vector>

using namespace tecplot::toolbox;
using std::vector;

void LoadData()
{
    Lock lockStart;
    /*
     * Create the Dataset (i.e data load step 2)
     */
    StringList varNames("X", "Y", NULL);
    TecUtilDataSetCreate("My DataSet", varNames.getRef(), TRUE);
    /*
     * Add a zone, and assign values to the variables. This amounts
     * to data load step 3.
     */

    FieldDataType_e dataTypes[] = {FieldDataType_Float,
                                   FieldDataType_Float
                                  };
    TecUtilDataSetAddZone("My Zone",
                          10, 10, 1, // zone dimensions
                          ZoneType_Ordered,
                          dataTypes);

    /*
     * Create temporary buffers to store the values in order to use
     * TecUtilDataValueArraySetByRef, rather than
     * TecUtilDataValueSetByRef. The former is preferred for
     * performance.
     */
    vector<float> xValues(100);
    vector<float> yValues(100);

    for (size_t j = 0; j < 10; j++)
        for (size_t i = 0; i < 10; i++)
        {
            size_t Node = j * 10 + i;
            xValues[Node] = static_cast<float>(i);
            yValues[Node] = static_cast<float>(j);
        }

    /*
     * Obtain writable references to variables 1 and 2 from zone 1,
     * and send the data to Tecplot.
     */
    FieldData_pa fdXVar = TecUtilDataValueGetWritableNativeRef(1, 1);
    FieldData_pa fdYVar = TecUtilDataValueGetWritableNativeRef(1, 2);

    TecUtilDataValueArraySetByRef(fdXVar,
                                  1, 100,
                                  &xValues[0]);
    TecUtilDataValueArraySetByRef(fdYVar,
                                  1, 100,
                                  &yValues[0]);

    /*
     * Finally inform the Tecplot engine that adding zones
     * is finished.
     */
    Set zoneSet(1);
    TecUtilStateChanged(StateChange_ZonesAdded,
                        (ArbParam_t)zoneSet.getRef());
    /*
     * This is not necessary to load data but if you want to see a
     * plot you need to switch the plottype to something reasonable.
    */
    TecUtilFrameSetPlotType(PlotType_Cartesian2D);
}
