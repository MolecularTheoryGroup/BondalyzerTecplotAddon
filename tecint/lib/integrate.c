#ifdef IRISX
#include <sys/types.h>
extern pid_t __vfork(void);
#endif

#include <math.h>
#include <limits.h>
#include <float.h>

#if !defined ADDON
#define ADDON
#endif /* ADDON */
#include "TECADDON.h"
#include "TASSERT.h"

#include "tecint.h"
#include "iresult.h"
#include "pctdone.h"
#include "zoneinfo.h"

static struct
{
    IntegrationResult_pa   IntegrationResult;
    VariableOption_e       VariableOption;
    Boolean_t              Axisymmetric;
    SymmetryVar_e          SymmetryVar;
    double                 SymmetryValue;
    EntIndex_t             ScalarVarNum;
    EntIndex_t             XVarNum;
    EntIndex_t             YVarNum;
    EntIndex_t             ZVarNum;
    IntegrateOver_e        IntegrateOver;
    Set_pa                 ZoneSet;
    IndexRange_t           IRange;
    IndexRange_t           JRange;
    IndexRange_t           KRange;
    Boolean_t              Absolute;
    Boolean_t              ExcludeBlanked;
} IntegrationValues;

static int ContainsNegativeVolumes;
static int ContainsPositiveVolumes;

static void GetCellNodes(LgIndex_t  CellNodes[8],
                         int       *NumCellNodes,
                         LgIndex_t  I,
                         LgIndex_t  J,
                         LgIndex_t  K)
/*
 * Get cell node indices.
 */
{
    LgIndex_t CellInd;

    REQUIRE(VALID_REF(CellNodes));
    REQUIRE(VALID_REF(NumCellNodes));
    REQUIRE(VALID_REF(CurZoneInfo));

    if (CurZoneInfo->ZoneType == ZoneType_Ordered)
    {
        *NumCellNodes = 2;
        CellNodes[0] = ((K - 1) * CurZoneInfo->JMax + (J - 1)) * CurZoneInfo->IMax + I;
        if (I < CurZoneInfo->IMax) /* Might be == IMax for face integrations. */
        {
            CellNodes[1] = CellNodes[0] + 1;
        }
        else /* Duplicate previous node (it shouldn't be used). */
        {
            CellNodes[1] = CellNodes[0];
        }
        if (J < CurZoneInfo->JMax)
        {
            CellNodes[2] = CellNodes[1] + CurZoneInfo->IMax;
            CellNodes[3] = CellNodes[0] + CurZoneInfo->IMax;
        }
        else
        {
            CellNodes[2] = CellNodes[1];
            CellNodes[3] = CellNodes[0];
        }
        if (K < CurZoneInfo->KMax)
        {
            CellNodes[4] = CellNodes[0] + CurZoneInfo->IMax * CurZoneInfo->JMax;
            CellNodes[5] = CellNodes[1] + CurZoneInfo->IMax * CurZoneInfo->JMax;
            CellNodes[6] = CellNodes[2] + CurZoneInfo->IMax * CurZoneInfo->JMax;
            CellNodes[7] = CellNodes[3] + CurZoneInfo->IMax * CurZoneInfo->JMax;
        }
        else
        {
            CellNodes[4] = CellNodes[0];
            CellNodes[5] = CellNodes[1];
            CellNodes[6] = CellNodes[2];
            CellNodes[7] = CellNodes[3];
        }

        if (CurZoneInfo->JMax > 1)
        {
            *NumCellNodes = 4;
        }

        if (CurZoneInfo->KMax > 1)
        {
            *NumCellNodes = 8;
        }
    }
    else
    {
        if (CurZoneInfo->ZoneType == ZoneType_FELineSeg)
        {
            *NumCellNodes = 2;
        }
        else if (CurZoneInfo->ZoneType == ZoneType_FETriangle)
        {
            *NumCellNodes = 3;
        }
        else if (CurZoneInfo->ZoneType == ZoneType_FEQuad)
        {
            *NumCellNodes = 4;
        }
        else if (CurZoneInfo->ZoneType == ZoneType_FETetra)
        {
            *NumCellNodes = 4;
        }
        else if (CurZoneInfo->ZoneType == ZoneType_FEBrick)
        {
            *NumCellNodes = 8;
        }
        for (CellInd = 0; CellInd < *NumCellNodes; CellInd++)
        {
            CellNodes[CellInd] = TecUtilDataNodeGetByRef(CurZoneInfo->NMap, J, CellInd + 1);
        }
    }
}


int CalculateCellContribution(double         *Numerator,
                              double         *Denominator,
                              FieldData_pa    Scalar,
                              FieldData_pa    Weight,
                              FieldData_pa    XComponent,
                              FieldData_pa    YComponent,
                              FieldData_pa    ZComponent,
                              IntegrateOver_e IntegrateOver,
                              LgIndex_t       I,
                              LgIndex_t       J,
                              LgIndex_t       K)
{
    /* Calculate the integral contribution(s) for a particular cell. For structured
       zones, this is cell (I,J,K) -> (I+1,J+1,K+1). For unstructured zones
       this is cell J, as described by NMap. */
    int FaceOffsets[4];
    int NumCellNodes;
    int NumUniqueCellNodes;
    int CellInd;
    int N;
    int TetraNodes[6][4] =
    /* For breaking a brick into six tetrahedra. */
    {
        {0, 2, 3, 7},
        {0, 4, 2, 7},
        {2, 7, 4, 6},
        {2, 4, 5, 6},
        {1, 5, 2, 4},
        {0, 1, 2, 4}
    };
    double Vector1[3];
    double Vector2[3];
    double CrossProduct[3];
    double Normal[3];
    double Tangential[3];
    double Length;
    double FaceArea1;
    double FaceArea2;
    double FaceArea;
    double Centroid1;
    double Centroid2;
    double Centroid;
    double Volume1;
    double Volume;
    double CellX[8];
    double CellY[8];
    double CellZ[8];
    double CellScalar[8];
    double CellWeight[8];
    double CellXComponent[8];
    double CellYComponent[8];
    double CellZComponent[8];
    double AveScalar;
    double AveWeight;
    double AveXComponent;
    double AveYComponent;
    double AveZComponent;
    LgIndex_t CellNodes[8];
    LgIndex_t UniqueCellNodes[8];
    Boolean_t IsBlanked = FALSE;
    Boolean_t IsOk = TRUE;

    REQUIRE(VALID_REF(Numerator));
    REQUIRE(VALID_REF(Denominator));
    REQUIRE(VALID_REF(Scalar) || Scalar == NULL);
    REQUIRE(VALID_REF(Weight) || Weight == NULL);
    REQUIRE(VALID_REF(XComponent) || XComponent == NULL);
    REQUIRE(VALID_REF(YComponent) || YComponent == NULL);
    REQUIRE(VALID_REF(ZComponent) || ZComponent == NULL);
    REQUIRE(CurZoneInfo->Zone > 0);

    *Numerator   = 0.0;
    *Denominator = 0.0;

    /* Get the nodes which make up the cell. */
    GetCellNodes(CellNodes, &NumCellNodes, I, J, K);

    /* Exclude edge/face/cell if blanked. */
    if (IntegrationValues.ExcludeBlanked)
    {
        if (CurZoneInfo->ZoneType == ZoneType_Ordered)
        {
            switch (IntegrateOver)
            {
                case IntegrateOver_Volume:
                    CellInd   = ((K - 1) * CurZoneInfo->JMax + J - 1) * CurZoneInfo->IMax + I;
                    IsBlanked = !TecUtilBlankingCheckIJKCell(CurZoneInfo->Zone, Planes_Volume, CellInd);
                    break;
                case IntegrateOver_IPlanes:
                    CellInd   = ((K - 1) * CurZoneInfo->JMax + J - 1) * CurZoneInfo->IMax + I;
                    IsBlanked = !TecUtilBlankingCheckIJKCell(CurZoneInfo->Zone, Planes_I, CellInd);
                    break;
                case IntegrateOver_JPlanes:
                    CellInd   = ((K - 1) * CurZoneInfo->JMax + J - 1) * CurZoneInfo->IMax + I;
                    IsBlanked = !TecUtilBlankingCheckIJKCell(CurZoneInfo->Zone, Planes_J, CellInd);
                    break;
                case IntegrateOver_KPlanes:
                    CellInd   = ((K - 1) * CurZoneInfo->JMax + J - 1) * CurZoneInfo->IMax + I;
                    IsBlanked = !TecUtilBlankingCheckIJKCell(CurZoneInfo->Zone, Planes_K, CellInd);
                    break;
                case IntegrateOver_ILines:
                    IsBlanked = !TecUtilBlankingCheckDataPoint(CurZoneInfo->Zone, CellNodes[0]) ||
                                !TecUtilBlankingCheckDataPoint(CurZoneInfo->Zone, CellNodes[1]);
                    break;
                case IntegrateOver_JLines:
                    IsBlanked = !TecUtilBlankingCheckDataPoint(CurZoneInfo->Zone, CellNodes[0]) ||
                                !TecUtilBlankingCheckDataPoint(CurZoneInfo->Zone, CellNodes[3]);
                    break;
                case IntegrateOver_KLines:
                    IsBlanked = !TecUtilBlankingCheckDataPoint(CurZoneInfo->Zone, CellNodes[0]) ||
                                !TecUtilBlankingCheckDataPoint(CurZoneInfo->Zone, CellNodes[4]);
                    break;
            }
        }
        else
        {
            IsBlanked = !TecUtilBlankingCheckFECell(CurZoneInfo->Zone, J);
        }
    }

    if (!IsBlanked)
    {
        Boolean_t IsScalarIntegrand;
        Boolean_t IsAveraged;

        /* Get integral form (scalar vs vector product, averaged vs not). */
        IsScalarIntegrand = (IntegrationValues.VariableOption == VariableOption_ScalarIntegral ||
                             IntegrationValues.VariableOption == VariableOption_Average ||
                             IntegrationValues.VariableOption == VariableOption_WeightedAverage);

        IsAveraged = (IntegrationValues.VariableOption == VariableOption_Average ||
                      IntegrationValues.VariableOption == VariableOption_WeightedAverage ||
                      IntegrationValues.VariableOption == VariableOption_VectorAverage);

        /* Get nodal values. */
        for (CellInd = 0; IsOk && CellInd < NumCellNodes; CellInd++)
        {
            CellXComponent[CellInd] = 0.0;
            CellYComponent[CellInd] = 0.0;
            CellZComponent[CellInd] = 0.0;

            if (CurZoneInfo->XVar)
            {
                CellX[CellInd] = TecUtilDataValueGetByRef(CurZoneInfo->XVar, CellNodes[CellInd]);
            }
            else
            {
                CellX[CellInd] = 0.0;
            }
            if (CurZoneInfo->YVar)
            {
                CellY[CellInd] = TecUtilDataValueGetByRef(CurZoneInfo->YVar, CellNodes[CellInd]);
            }
            else
            {
                CellY[CellInd] = 0.0;
            }
            if (CurZoneInfo->ZVar)
            {
                CellZ[CellInd] = TecUtilDataValueGetByRef(CurZoneInfo->ZVar, CellNodes[CellInd]);
            }
            else
            {
                CellZ[CellInd] = 0.0;
            }
            if (Scalar)
            {
                CellScalar[CellInd] = TecUtilDataValueGetByRef(Scalar, CellNodes[CellInd]);
            }
            else
            {
                CellScalar[CellInd] = 1.0;
            }
            if (Weight)
            {
                CellWeight[CellInd] = TecUtilDataValueGetByRef(Weight, CellNodes[CellInd]);
            }
            else
            {
                CellWeight[CellInd] = 1.0;
            }
            if (XComponent)
            {
                CellXComponent[CellInd] = TecUtilDataValueGetByRef(XComponent, CellNodes[CellInd]);
            }
            if (YComponent)
            {
                CellYComponent[CellInd] = TecUtilDataValueGetByRef(YComponent, CellNodes[CellInd]);
            }
            if (ZComponent)
            {
                CellZComponent[CellInd] = TecUtilDataValueGetByRef(ZComponent, CellNodes[CellInd]);
            }
        }

        if (IsOk)
        {
            /* For triangular cells, set node 4 equal to node 3 */
            if (CurZoneInfo->ZoneType == ZoneType_FETriangle)
            {
                CellNodes     [3] = CellNodes     [2];
                CellX         [3] = CellX         [2];
                CellY         [3] = CellY         [2];
                CellZ         [3] = CellZ         [2];
                CellScalar    [3] = CellScalar    [2];
                CellWeight    [3] = CellWeight    [2];
                CellXComponent[3] = CellXComponent[2];
                CellYComponent[3] = CellYComponent[2];
                CellZComponent[3] = CellZComponent[2];
            }

            switch (IntegrateOver)
            {
                case IntegrateOver_Volume:
                    /* Calculate the cell volume. */
                    if (CurZoneInfo->ZoneType == ZoneType_FETetra)
                    {
                        Vector1[0] = CellX[1] - CellX[0];
                        Vector1[1] = CellY[1] - CellY[0];
                        Vector1[2] = CellZ[1] - CellZ[0];
                        Vector2[0] = CellX[2] - CellX[0];
                        Vector2[1] = CellY[2] - CellY[0];
                        Vector2[2] = CellZ[2] - CellZ[0];
                        CrossProduct[0] = Vector1[1] * Vector2[2] -
                                          Vector2[1] * Vector1[2];
                        CrossProduct[1] = Vector1[2] * Vector2[0] -
                                          Vector2[2] * Vector1[0];
                        CrossProduct[2] = Vector1[0] * Vector2[1] -
                                          Vector2[0] * Vector1[1];
                        Vector1[0] = CellX[3] - CellX[0];
                        Vector1[1] = CellY[3] - CellY[0];
                        Vector1[2] = CellZ[3] - CellZ[0];
                        Volume = (CrossProduct[0] * Vector1[0] +
                                  CrossProduct[1] * Vector1[1] +
                                  CrossProduct[2] * Vector1[2]) / 6.0;
                    }
                    else
                    {
                        /* Add the volumes of six tetrahedra within the brick. */
                        int n;

                        Volume = 0.0;
                        for (n = 0; n < 6; n++)
                        {
                            Vector1[0] = CellX[TetraNodes[n][1]] - CellX[TetraNodes[n][0]];
                            Vector1[1] = CellY[TetraNodes[n][1]] - CellY[TetraNodes[n][0]];
                            Vector1[2] = CellZ[TetraNodes[n][1]] - CellZ[TetraNodes[n][0]];
                            Vector2[0] = CellX[TetraNodes[n][2]] - CellX[TetraNodes[n][0]];
                            Vector2[1] = CellY[TetraNodes[n][2]] - CellY[TetraNodes[n][0]];
                            Vector2[2] = CellZ[TetraNodes[n][2]] - CellZ[TetraNodes[n][0]];
                            CrossProduct[0] = Vector1[1] * Vector2[2] -
                                              Vector2[1] * Vector1[2];
                            CrossProduct[1] = Vector1[2] * Vector2[0] -
                                              Vector2[2] * Vector1[0];
                            CrossProduct[2] = Vector1[0] * Vector2[1] -
                                              Vector2[0] * Vector1[1];
                            Vector1[0] = CellX[TetraNodes[n][3]] - CellX[TetraNodes[n][0]];
                            Vector1[1] = CellY[TetraNodes[n][3]] - CellY[TetraNodes[n][0]];
                            Vector1[2] = CellZ[TetraNodes[n][3]] - CellZ[TetraNodes[n][0]];
                            Volume1 = CrossProduct[0] * Vector1[0] +
                                      CrossProduct[1] * Vector1[1] +
                                      CrossProduct[2] * Vector1[2];
                            Volume += Volume1;
                        }
                        Volume /= 6.0;
                    }
                    if (Volume < 0.0)
                    {
                        if (IntegrationValues.Absolute)
                        {
                            Volume = -Volume;
                        }
                        else
                        {
                            ContainsNegativeVolumes = 1;
                        }
                    }
                    else if (Volume > 0.0)
                    {
                        ContainsPositiveVolumes = 1;
                    }
                    if (IntegrationValues.VariableOption == VariableOption_LengthAreaVolume)
                    {
                        *Numerator = Volume;
                    }
                    else
                    {
                        /* Get the average scalar and weight. */
                        NumUniqueCellNodes = 0;
                        AveScalar = 0.0;
                        AveWeight = 0.0;

                        for (CellInd = 0; CellInd < NumCellNodes; CellInd++)
                        {
                            for (N = 0; N < NumUniqueCellNodes; N++)
                            {
                                if (CellNodes[CellInd] == UniqueCellNodes[N])
                                {
                                    break;
                                }
                            }
                            if (N == NumUniqueCellNodes) /* Found a unique node. */
                            {
                                AveScalar += CellScalar[CellInd];
                                AveWeight += CellWeight[CellInd];
                                UniqueCellNodes[NumUniqueCellNodes++] = CellNodes[CellInd];
                            }
                        }

                        AveScalar /= (double)NumUniqueCellNodes;
                        AveWeight /= (double)NumUniqueCellNodes;

                        /* The cell contribution is the volume multiplied by the average cell value. */
                        *Numerator = Volume * AveScalar * AveWeight;

                        if (IsAveraged)
                        {
                            *Denominator = Volume * AveWeight;
                        }
                    }
                    break;
                case IntegrateOver_IPlanes: /* These can be scalar or vector dotted with face normal. */
                case IntegrateOver_JPlanes:
                case IntegrateOver_KPlanes:
                    if (IntegrateOver == IntegrateOver_IPlanes)
                    {
                        FaceOffsets[0] = 0;
                        FaceOffsets[1] = 3;
                        FaceOffsets[2] = 7;
                        FaceOffsets[3] = 4;
                    }
                    else if (IntegrateOver == IntegrateOver_JPlanes)
                    {
                        FaceOffsets[0] = 0;
                        FaceOffsets[1] = 4;
                        FaceOffsets[2] = 5;
                        FaceOffsets[3] = 1;
                    }
                    else
                    {
                        FaceOffsets[0] = 0;
                        FaceOffsets[1] = 1;
                        FaceOffsets[2] = 2;
                        FaceOffsets[3] = 3;
                    }

                    /* The face area is calculated with cross products of adjacent edge vectors. */
                    Vector1[0] = CellX[FaceOffsets[1]] - CellX[FaceOffsets[0]];
                    Vector1[1] = CellY[FaceOffsets[1]] - CellY[FaceOffsets[0]];
                    Vector1[2] = CellZ[FaceOffsets[1]] - CellZ[FaceOffsets[0]];
                    Vector2[0] = CellX[FaceOffsets[3]] - CellX[FaceOffsets[0]];
                    Vector2[1] = CellY[FaceOffsets[3]] - CellY[FaceOffsets[0]];
                    Vector2[2] = CellZ[FaceOffsets[3]] - CellZ[FaceOffsets[0]];
                    CrossProduct[0] = Vector1[1] * Vector2[2] -
                                      Vector2[1] * Vector1[2];
                    CrossProduct[1] = Vector1[2] * Vector2[0] -
                                      Vector2[2] * Vector1[0];
                    CrossProduct[2] = Vector1[0] * Vector2[1] -
                                      Vector2[0] * Vector1[1];
                    FaceArea1 = sqrt(CrossProduct[0] * CrossProduct[0] +
                                     CrossProduct[1] * CrossProduct[1] +
                                     CrossProduct[2] * CrossProduct[2]);
                    Normal[0] = CrossProduct[0];
                    Normal[1] = CrossProduct[1];
                    Normal[2] = CrossProduct[2];
                    Vector1[0] = CellX[FaceOffsets[3]] - CellX[FaceOffsets[2]];
                    Vector1[1] = CellY[FaceOffsets[3]] - CellY[FaceOffsets[2]];
                    Vector1[2] = CellZ[FaceOffsets[3]] - CellZ[FaceOffsets[2]];
                    Vector2[0] = CellX[FaceOffsets[1]] - CellX[FaceOffsets[2]];
                    Vector2[1] = CellY[FaceOffsets[1]] - CellY[FaceOffsets[2]];
                    Vector2[2] = CellZ[FaceOffsets[1]] - CellZ[FaceOffsets[2]];
                    CrossProduct[0] = Vector1[1] * Vector2[2] -
                                      Vector2[1] * Vector1[2];
                    CrossProduct[1] = Vector1[2] * Vector2[0] -
                                      Vector2[2] * Vector1[0];
                    CrossProduct[2] = Vector1[0] * Vector2[1] -
                                      Vector2[0] * Vector1[1];
                    FaceArea2 = sqrt(CrossProduct[0] * CrossProduct[0] +
                                     CrossProduct[1] * CrossProduct[1] +
                                     CrossProduct[2] * CrossProduct[2]);
                    FaceArea = 0.5 * (FaceArea1 + FaceArea2);

                    if (IntegrationValues.Axisymmetric)
                    {
                        /* Calculate the face centroid. */
                        if (IntegrationValues.SymmetryVar == SymmetryVar_X)
                        {
                            /* Calculate the Distance of the centroid from the given X value. */
                            Centroid1 = (CellX[FaceOffsets[0]] + CellX[FaceOffsets[1]] +
                                         CellX[FaceOffsets[3]]) / 3.0 - IntegrationValues.SymmetryValue;
                            Centroid2 = (CellX[FaceOffsets[1]] + CellX[FaceOffsets[2]] +
                                         CellX[FaceOffsets[3]]) / 3.0 - IntegrationValues.SymmetryValue;
                        }
                        else
                        {
                            /* Calculate the Distance of the centroid from the given X value. */
                            Centroid1 = (CellY[FaceOffsets[0]] + CellY[FaceOffsets[1]] +
                                         CellY[FaceOffsets[3]]) / 3.0 - IntegrationValues.SymmetryValue;
                            Centroid2 = (CellY[FaceOffsets[1]] + CellY[FaceOffsets[2]] +
                                         CellY[FaceOffsets[3]]) / 3.0 - IntegrationValues.SymmetryValue;
                        }

                        /* Now area-weight the two centroids, and multiply by two pi. */
                        if (fabs(FaceArea) > FLT_MIN)
                        {
                            Centroid = PI * (Centroid1 * FaceArea1 + Centroid2 * FaceArea2) / FaceArea;
                        }
                        else
                        {
                            Centroid = 0.0;
                        }
                    }
                    else
                    {
                        Centroid = 1.0;
                    }

                    /* Multiply the face area by 2 * PI * Centroid (or unity if not axisymmetric). */
                    FaceArea *= Centroid;

                    if (FaceArea < 0.0)
                    {
                        if (IntegrationValues.Absolute)
                        {
                            FaceArea = -FaceArea;
                        }
                        else
                        {
                            ContainsNegativeVolumes = 1;
                        }
                    }
                    else if (FaceArea > 0.0)
                    {
                        ContainsPositiveVolumes = 1;
                    }
                    Normal[0] += CrossProduct[0];
                    Normal[1] += CrossProduct[1];
                    Normal[2] += CrossProduct[2];
                    if (IntegrationValues.VariableOption == VariableOption_LengthAreaVolume)
                    {
                        *Numerator = FaceArea;
                    }
                    else if (IsScalarIntegrand)
                    {
                        /* Get the average scalar and weight. */
                        NumUniqueCellNodes = 0;
                        AveScalar = 0.0;
                        AveWeight = 0.0;

                        for (CellInd = 0; CellInd < 4; CellInd++)
                        {
                            for (N = 0; N < NumUniqueCellNodes; N++)
                            {
                                if (CellNodes[FaceOffsets[CellInd]] == UniqueCellNodes[N])
                                {
                                    break;
                                }
                            }

                            if (N == NumUniqueCellNodes) /* Found a unique node. */
                            {
                                AveScalar += CellScalar[FaceOffsets[CellInd]];
                                AveWeight += CellWeight[FaceOffsets[CellInd]];
                                UniqueCellNodes[NumUniqueCellNodes++] = CellNodes[FaceOffsets[CellInd]];
                            }
                        }

                        AveScalar /= (double)NumUniqueCellNodes;
                        AveWeight /= (double)NumUniqueCellNodes;

                        /* The face contribution is the area multiplied by the average face value. */
                        *Numerator = FaceArea * AveScalar * AveWeight;

                        if (IsAveraged)
                        {
                            *Denominator = FaceArea * AveWeight;
                        }
                    }
                    else /* Vector-dot-normal integration, including Forces & Moments. */
                    {
                        double DotProduct;

                        /* normalise the average face normal. A normalized normal -- cool. */
                        Length = sqrt(Normal[0] * Normal[0] +
                                      Normal[1] * Normal[1] +
                                      Normal[2] * Normal[2]);
                        if (Length > FLT_MIN)
                        {
                            Normal[0] /= Length;
                            Normal[1] /= Length;
                            Normal[2] /= Length;
                        }
                        /* Calculate the average face value of each component from the unique nodes only,
                         * ignoring repeated nodes. */
                        NumUniqueCellNodes = 0;
                        AveXComponent = 0.0;
                        AveYComponent = 0.0;
                        AveZComponent = 0.0;
                        AveScalar     = 0.0;
                        for (CellInd = 0; CellInd < 4; CellInd++)
                        {
                            for (N = 0; N < NumUniqueCellNodes; N++)
                            {
                                if (CellNodes[FaceOffsets[CellInd]] == UniqueCellNodes[N])
                                {
                                    break;
                                }
                            }
                            if (N == NumUniqueCellNodes)
                            {
                                /* Found a unique node. */
                                AveXComponent += CellXComponent[FaceOffsets[CellInd]];
                                AveYComponent += CellYComponent[FaceOffsets[CellInd]];
                                AveZComponent += CellZComponent[FaceOffsets[CellInd]];
                                AveScalar     += CellScalar[FaceOffsets[CellInd]];
                                UniqueCellNodes[NumUniqueCellNodes++] = CellNodes[FaceOffsets[CellInd]];
                            }
                        }
                        AveXComponent /= (double)NumUniqueCellNodes;
                        AveYComponent /= (double)NumUniqueCellNodes;
                        AveZComponent /= (double)NumUniqueCellNodes;
                        AveScalar     /= (double)NumUniqueCellNodes;

                        /* General vector-dot-normal integration. */
                        DotProduct = (AveXComponent * Normal[0] +
                                      AveYComponent * Normal[1] +
                                      AveZComponent * Normal[2]);
                        *Numerator = FaceArea * AveScalar * DotProduct;
                        if (IntegrationValues.VariableOption == VariableOption_VectorAverage)
                        {
                            *Denominator = FaceArea;
                        }
                    }
                    break;
                case IntegrateOver_ILines: /* These can be scalar integrations or vector dotted with a normal or the tangential. */
                case IntegrateOver_JLines: /* The calculation of the normal is done so as to be consistent with 3D face normals, */
                case IntegrateOver_KLines: /* which point in the plus-index direction. Accordingly, left-normals are used for    */
                    /* I lines, right-normals are used for J-lines, and K-line normals are set to zero.   */
                    if (IntegrateOver == IntegrateOver_ILines)
                    {
                        CellInd = 1;
                    }
                    else if (IntegrateOver == IntegrateOver_JLines)
                    {
                        CellInd = 3;
                    }
                    else
                    {
                        CellInd = 4;
                    }

                    Tangential[0] = CellX[CellInd] - CellX[0];
                    Tangential[1] = CellY[CellInd] - CellY[0];
                    Tangential[2] = CellZ[CellInd] - CellZ[0];

                    AveScalar     = 0.5 * (CellScalar[0] + CellScalar[CellInd]);
                    AveWeight     = 0.5 * (CellWeight[0] + CellWeight[CellInd]);
                    AveXComponent = 0.5 * (CellXComponent[0] + CellXComponent[CellInd]);
                    AveYComponent = 0.5 * (CellYComponent[0] + CellYComponent[CellInd]);
                    AveZComponent = 0.5 * (CellZComponent[0] + CellZComponent[CellInd]);

                    if (IntegrationValues.Axisymmetric)
                    {
                        if (IntegrationValues.SymmetryVar == SymmetryVar_X)
                        {
                            Centroid = PI * (CellX[0] + CellX[CellInd] - 2.0 * IntegrationValues.SymmetryValue);
                        }
                        else
                        {
                            Centroid = PI * (CellY[0] + CellY[CellInd] - 2.0 * IntegrationValues.SymmetryValue);
                        }
                    }
                    else
                    {
                        Centroid = 1.0;
                    }

                    Length = Centroid * sqrt(Tangential[0] * Tangential[0] +
                                             Tangential[1] * Tangential[1] +
                                             Tangential[2] * Tangential[2]);

                    if (IntegrationValues.VariableOption == VariableOption_LengthAreaVolume)
                    {
                        *Numerator = Length;
                    }
                    else if (IsScalarIntegrand)
                    {
                        *Numerator = Length * AveWeight * AveScalar;

                        if (IsAveraged)
                        {
                            *Denominator = Length * AveWeight;
                        }
                    }
                    else
                    {
                        /* Vector-type integrands. */
                        if (IntegrationValues.VariableOption == VariableOption_VectorDotTangential)
                        {
                            *Numerator = Tangential[0] * AveXComponent +
                                         Tangential[1] * AveYComponent +
                                         Tangential[2] * AveZComponent;
                        }
                        else
                        {
                            /* Vector-dot-normal integrations. */
                            double NormalLength;
                            double DotProduct;

                            if (IntegrateOver == IntegrateOver_ILines)
                            {
                                Normal[0] = -Tangential[1];
                                Normal[1] = Tangential[0];
                            }
                            else if (IntegrateOver == IntegrateOver_JLines)
                            {
                                Normal[0] = Tangential[1];
                                Normal[1] = -Tangential[0];
                            }
                            else
                            {
                                Normal[0] = 0.0;
                                Normal[1] = 0.0;
                            }
                            Normal[2]    = 0.0;
                            NormalLength = sqrt(Normal[0] * Normal[0] +
                                                Normal[1] * Normal[1]);
                            if (NormalLength > FLT_MIN)
                            {
                                Normal[0] /= NormalLength;
                                Normal[1] /= NormalLength;
                            }

                            /* General vector-dot-normal integral */
                            DotProduct = (Normal[0] * AveXComponent +
                                          Normal[1] * AveYComponent +
                                          Normal[2] * AveZComponent);

                            *Numerator = Length * AveScalar * DotProduct;
                            if (IntegrationValues.VariableOption == VariableOption_VectorAverage)
                            {
                                *Denominator = Length;
                            }
                        }
                    }
                    break;
                default:
                    IsOk = FALSE;
            }
        } /* IsOk */
    } /* !IsBlanked */

    return(IsOk);
}

ReturnStatus_e PerformSingleIntegration(double          *Numerator,
                                        double          *Denominator,
                                        FieldData_pa     Scalar,
                                        FieldData_pa     Weight,
                                        FieldData_pa     XComponent,
                                        FieldData_pa     YComponent,
                                        FieldData_pa     ZComponent,
                                        IntegrateOver_e  IntegrateOver,
                                        int              IMin,
                                        int              JMin,
                                        int              KMin,
                                        int              IMax,
                                        int              JMax,
                                        int              KMax,
                                        int              PercentDoneMinimum,
                                        int              PercentDoneMaximum)
{
    int            I;
    int            J;
    int            K;
    int            PercentDone;
    int            DeltaPercentDone = PercentDoneMaximum - PercentDoneMinimum;
    double         NumContrib;
    double         DenContrib;
    ReturnStatus_e ReturnStatus = ReturnStatus_OK;

    *Numerator   = 0.0;
    *Denominator = 0.0;

    switch (IntegrateOver)
    {
        case IntegrateOver_Volume: /* Must be 3D, can be structured or unstructured */
            if (CurZoneInfo->ZoneType == ZoneType_Ordered)
            {
                for (K = KMin; ReturnStatus == ReturnStatus_OK && K < KMax; K++)
                {
                    for (J = JMin; ReturnStatus == ReturnStatus_OK && J < JMax; J++)
                    {
                        PercentDone = PercentDoneMinimum + DeltaPercentDone *
                                      ((K - KMin) * (JMax - JMin) + J - JMin) /
                                      ((KMax - KMin) * (JMax - JMin));
                        if (!CheckPercentDone(PercentDone))
                        {
                            ReturnStatus = ReturnStatus_Canceled;
                        }
                        for (I = IMin; ReturnStatus == ReturnStatus_OK && I < IMax; I++)
                        {
                            if (CalculateCellContribution(&NumContrib,
                                                          &DenContrib,
                                                          Scalar,
                                                          Weight,
                                                          XComponent,
                                                          YComponent,
                                                          ZComponent,
                                                          IntegrateOver,
                                                          I,
                                                          J,
                                                          K))
                            {
                                *Numerator   += NumContrib;
                                *Denominator += DenContrib;
                            }
                            else
                            {
                                ReturnStatus = ReturnStatus_Error;
                            }
                        } /* I */
                    } /* J */
                } /* K */
            } /* Ordered zone. */
            else if (CurZoneInfo->ZoneType == ZoneType_FEBrick ||
                     CurZoneInfo->ZoneType == ZoneType_FETetra)
            {
                for (J = JMin; ReturnStatus == ReturnStatus_OK && J <= JMax; J++)
                {
                    if (J % 50 == 0)
                    {
                        PercentDone = PercentDoneMinimum + DeltaPercentDone *
                                      (J - JMin) / (JMax - JMin);
                        if (!CheckPercentDone(PercentDone))
                        {
                            ReturnStatus = ReturnStatus_Canceled;
                        }
                    }
                    if (CalculateCellContribution(&NumContrib,
                                                  &DenContrib,
                                                  Scalar,
                                                  Weight,
                                                  XComponent,
                                                  YComponent,
                                                  ZComponent,
                                                  IntegrateOver,
                                                  IMin,
                                                  J,
                                                  KMin))
                    {
                        *Numerator   += NumContrib;
                        *Denominator += DenContrib;
                    }
                    else
                    {
                        ReturnStatus = ReturnStatus_Error;
                    }
                    /* Cannot be forces & moments, so only one result */
                }
            } /* FE Zone. */
            break;
        case IntegrateOver_IPlanes: /* Must be 3D structured zone */
            for (K = KMin; ReturnStatus == ReturnStatus_OK && K < KMax; K++)
            {
                PercentDone = PercentDoneMinimum + DeltaPercentDone *
                              (K - KMin) / (KMax - KMin);
                if (!CheckPercentDone(PercentDone))
                {
                    ReturnStatus = ReturnStatus_Canceled;
                }

                for (J = JMin; ReturnStatus == ReturnStatus_OK && J < JMax; J++)
                {
                    if (CalculateCellContribution(&NumContrib,
                                                  &DenContrib,
                                                  Scalar,
                                                  Weight,
                                                  XComponent,
                                                  YComponent,
                                                  ZComponent,
                                                  IntegrateOver,
                                                  IMin,
                                                  J,
                                                  K))
                    {
                        *Numerator   += NumContrib;
                        *Denominator += DenContrib;
                    }
                    else
                    {
                        ReturnStatus = ReturnStatus_Error;
                    }
                }
            }
            break;
        case IntegrateOver_JPlanes: /* Must be 3D structured zone */
            for (K = KMin; ReturnStatus == ReturnStatus_OK && K < KMax; K++)
            {
                PercentDone = PercentDoneMinimum + DeltaPercentDone *
                              (K - KMin) / (KMax - KMin);
                if (!CheckPercentDone(PercentDone))
                {
                    ReturnStatus = ReturnStatus_Canceled;
                }
                for (I = IMin; ReturnStatus == ReturnStatus_OK && I < IMax; I++)
                {
                    if (CalculateCellContribution(&NumContrib,
                                                  &DenContrib,
                                                  Scalar,
                                                  Weight,
                                                  XComponent,
                                                  YComponent,
                                                  ZComponent,
                                                  IntegrateOver,
                                                  I,
                                                  JMin,
                                                  K))
                    {
                        *Numerator   += NumContrib;
                        *Denominator += DenContrib;
                    }
                    else
                    {
                        ReturnStatus = ReturnStatus_Error;
                    }
                }
            }
            break;
        case IntegrateOver_KPlanes: /* Can be 2D, can be unstructured */
            if (CurZoneInfo->ZoneType == ZoneType_Ordered)
            {
                for (J = JMin; ReturnStatus == ReturnStatus_OK && J < JMax; J++)
                {
                    PercentDone = PercentDoneMinimum + DeltaPercentDone *
                                  (J - JMin) / (JMax - JMin);
                    if (!CheckPercentDone(PercentDone))
                    {
                        ReturnStatus = ReturnStatus_Canceled;
                    }

                    for (I = IMin; ReturnStatus == ReturnStatus_OK && I < IMax; I++)
                    {
                        if (CalculateCellContribution(&NumContrib,
                                                      &DenContrib,
                                                      Scalar,
                                                      Weight,
                                                      XComponent,
                                                      YComponent,
                                                      ZComponent,
                                                      IntegrateOver,
                                                      I,
                                                      J,
                                                      KMin))
                        {
                            *Numerator   += NumContrib;
                            *Denominator += DenContrib;
                        }
                        else
                        {
                            ReturnStatus = ReturnStatus_Error;
                        }
                    }
                }
            }
            else
            {
                for (J = JMin; ReturnStatus == ReturnStatus_OK && J <= JMax; J++)
                {
                    if (J % 50 == 0)
                    {
                        PercentDone = PercentDoneMinimum + DeltaPercentDone *
                                      (J - JMin) / (JMax - JMin);
                        if (!CheckPercentDone(PercentDone))
                        {
                            ReturnStatus = ReturnStatus_Canceled;
                        }
                    }

                    if (CalculateCellContribution(&NumContrib,
                                                  &DenContrib,
                                                  Scalar,
                                                  Weight,
                                                  XComponent,
                                                  YComponent,
                                                  ZComponent,
                                                  IntegrateOver,
                                                  IMin,
                                                  J,
                                                  KMin))
                    {
                        *Numerator   += NumContrib;
                        *Denominator += DenContrib;
                    }
                    else
                    {
                        ReturnStatus = ReturnStatus_Error;
                    }
                }
            }
            break;
        case IntegrateOver_ILines: /* Must be FELineSeg or structured 1D, 2D, or 3D. */
            if (CurZoneInfo->ZoneType == ZoneType_FELineSeg)
            {
                for (J = JMin; ReturnStatus == ReturnStatus_OK && J <= JMax; J++)
                {
                    if (J % 50 == 0)
                    {
                        PercentDone = PercentDoneMinimum + DeltaPercentDone *
                                      (J - JMin) / (JMax - JMin);
                        if (!CheckPercentDone(PercentDone))
                        {
                            ReturnStatus = ReturnStatus_Canceled;
                        }
                    }

                    if (CalculateCellContribution(&NumContrib,
                                                  &DenContrib,
                                                  Scalar,
                                                  Weight,
                                                  XComponent,
                                                  YComponent,
                                                  ZComponent,
                                                  IntegrateOver,
                                                  IMin,
                                                  J,
                                                  KMin))
                    {
                        *Numerator   += NumContrib;
                        *Denominator += DenContrib;
                    }
                    else
                    {
                        ReturnStatus = ReturnStatus_Error;
                    }
                }
            }
            else
            {
                for (I = IMin; ReturnStatus == ReturnStatus_OK && I < IMax; I++)
                {
                    if (I % 50 == 0)
                    {
                        PercentDone = PercentDoneMinimum + DeltaPercentDone *
                                      (I - IMin) / (IMax - IMin);
                        if (!CheckPercentDone(PercentDone))
                        {
                            ReturnStatus = ReturnStatus_Canceled;
                        }
                    }

                    if (CalculateCellContribution(&NumContrib,
                                                  &DenContrib,
                                                  Scalar,
                                                  Weight,
                                                  XComponent,
                                                  YComponent,
                                                  ZComponent,
                                                  IntegrateOver,
                                                  I,
                                                  JMin,
                                                  KMin))
                    {
                        *Numerator   += NumContrib;
                        *Denominator += DenContrib;
                    }
                    else
                    {
                        ReturnStatus = ReturnStatus_Error;
                    }
                }
            }
            break;
        case IntegrateOver_JLines: /* Must be 2D or 3D structured. */
            for (J = JMin; ReturnStatus == ReturnStatus_OK && J < JMax; J++)
            {
                if (J % 50 == 0)
                {
                    PercentDone = PercentDoneMinimum + DeltaPercentDone *
                                  (J - JMin) / (JMax - JMin);
                    if (!CheckPercentDone(PercentDone))
                    {
                        ReturnStatus = ReturnStatus_Canceled;
                    }
                }

                if (CalculateCellContribution(&NumContrib,
                                              &DenContrib,
                                              Scalar,
                                              Weight,
                                              XComponent,
                                              YComponent,
                                              ZComponent,
                                              IntegrateOver,
                                              IMin,
                                              J,
                                              KMin))
                {
                    *Numerator   += NumContrib;
                    *Denominator += DenContrib;
                }
                else
                {
                    ReturnStatus = ReturnStatus_Error;
                }
            }
            break;
        case IntegrateOver_KLines: /* Must be 3D structured. */
            for (K = KMin; ReturnStatus == ReturnStatus_OK && K < KMax; K++)
            {
                if (K % 50 == 0)
                {
                    PercentDone = PercentDoneMinimum + DeltaPercentDone *
                                  (K - KMin) / (KMax - KMin);
                    if (!CheckPercentDone(PercentDone))
                    {
                        ReturnStatus = ReturnStatus_Canceled;
                    }
                }

                if (CalculateCellContribution(&NumContrib,
                                              &DenContrib,
                                              Scalar,
                                              Weight,
                                              XComponent,
                                              YComponent,
                                              ZComponent,
                                              IntegrateOver,
                                              IMin,
                                              JMin,
                                              K))
                {
                    *Numerator   += NumContrib;
                    *Denominator += DenContrib;
                }
                else
                {
                    ReturnStatus = ReturnStatus_Error;
                }
            }
            break;
        default:
            TecUtilDialogErrMsg("PerformSingleIntegration: Unrecognized\n"
                                "integration option (internal error).");
            ReturnStatus = ReturnStatus_Error;
            break;
    }

    ENSURE(VALID_ENUM(ReturnStatus, ReturnStatus_e));
    return(ReturnStatus);
}

ReturnStatus_e DoIntegration(void)
{
    /* Perform an integration over the zones set in IntegrationValues.ZoneSet, and store the results to
       Numerator and Denominator, which must already have been allocated as an array of pointers
       (should eventually replace this with some self-describing structure allocated here).
       These should be sized as
          double *Numerator[MaxZoneNum], and each element of Numerator and Denominator should be
          double Numerator[MaxZoneNum][NumberOfValuesPerZone]. See the below switch(IntegrateOver) statement
       for hints as to how to calculate NumberOfValuesPerZone. It isn't pretty.

       Face normals are taken in the + index directions (e.g. normals to I=constant planes point in the
       +I direction). For consistency, line left normals are used.
    */

    char            PercentDialogText[256];
    int             I;
    int             J;
    int             K;
    int             NumberOfIntegrations;
    int             IntegrationNumber;
    int             PercentDoneMinimum;
    int             PercentDoneMaximum;
    double          Numerator;
    double          Denominator;
    double          NumeratorTotal;
    double          DenominatorTotal;
    EntIndex_t      MaxZoneNum;
    EntIndex_t      ZoneNum;
    LgIndex_t       IMin;
    LgIndex_t       JMin;
    LgIndex_t       KMin;
    LgIndex_t       IMax;
    LgIndex_t       JMax;
    LgIndex_t       KMax;
    FieldData_pa    ScalarVar;
    FieldData_pa    WeightingVar;
    FieldData_pa    XComponent;
    FieldData_pa    YComponent;
    FieldData_pa    ZComponent;
    FrameMode_e     FrameMode = TecUtilFrameGetMode();
    Boolean_t       UsesScalarVariable;
    IntegrateOver_e IntegrateOver;
    ReturnStatus_e  ReturnStatus = ReturnStatus_OK;

    REQUIRE(FrameMode == Frame_XY || FrameMode == Frame_TwoD || FrameMode == Frame_ThreeD);

    UsesScalarVariable = !(IntegrationValues.VariableOption == VariableOption_LengthAreaVolume ||
                           IntegrationValues.VariableOption == VariableOption_VectorDotNormal ||
                           IntegrationValues.VariableOption == VariableOption_VectorAverage ||
                           IntegrationValues.VariableOption == VariableOption_VectorDotTangential);

    REQUIRE(IMPLICATION(UsesScalarVariable, TecUtilVarIsEnabled(IntegrationValues.ScalarVarNum)));

    /* Track whether this integration is over a left-handed or right-handed
    domain, or (gasp!) both. */
    ContainsNegativeVolumes = 0;
    ContainsPositiveVolumes = 0;

    /* Zero out the overall total(s). */
    TecUtilDataSetGetInfo(NULL, &MaxZoneNum, NULL);
    NumeratorTotal = 0.0;
    DenominatorTotal = 0.0;

    FrameMode = TecUtilFrameGetMode();

    /* Loop over all assigned zones. */
    ZoneNum = (EntIndex_t)TecUtilSetGetNextMember(IntegrationValues.ZoneSet,
                                                  TECUTILSETNOTMEMBER);
    while (ReturnStatus == ReturnStatus_OK && ZoneNum != (EntIndex_t)TECUTILSETNOTMEMBER)
    {
        double ZoneNumeratorTotal   = 0.0;
        double ZoneDenominatorTotal = 0.0;

        if (TecUtilZoneIsFiniteElement(ZoneNum) &&
            TecUtilZoneGetType(ZoneNum) != ZoneType_FELineSeg &&
            FrameMode != Frame_TwoD &&
            FrameMode != Frame_ThreeD)
        {
            TecUtilDialogErrMsg("Please set plot type to 2D or 3D Cartesian\n"
                                "to integrate finite-element area/volume zones.");
            ReturnStatus = ReturnStatus_Error;
        }

        if (ReturnStatus == ReturnStatus_OK)
        {
            SetCurZoneInfo(ZoneNum);

            if (IntegrationValues.VariableOption == VariableOption_VectorDotNormal ||
                IntegrationValues.VariableOption == VariableOption_VectorAverage ||
                IntegrationValues.VariableOption == VariableOption_VectorDotTangential)
            {
                if (IntegrationValues.XVarNum &&
                    TecUtilVarIsEnabled(IntegrationValues.XVarNum))
                {
                    XComponent = TecUtilDataValueGetRef(ZoneNum, IntegrationValues.XVarNum);
                }
                else
                {
                    XComponent = NULL;
                }
                if (IntegrationValues.YVarNum &&
                    TecUtilVarIsEnabled(IntegrationValues.YVarNum))
                {
                    YComponent = TecUtilDataValueGetRef(ZoneNum, IntegrationValues.YVarNum);
                }
                else
                {
                    YComponent = NULL;
                }
                if (IntegrationValues.ZVarNum &&
                    TecUtilVarIsEnabled(IntegrationValues.ZVarNum))
                {
                    ZComponent = TecUtilDataValueGetRef(ZoneNum, IntegrationValues.ZVarNum);
                }
                else
                {
                    ZComponent = NULL;
                }
                WeightingVar = NULL;
            }
            else if (IntegrationValues.VariableOption == VariableOption_WeightedAverage)
            {
                XComponent = NULL;
                YComponent = NULL;
                ZComponent = NULL;

                if (IntegrationValues.XVarNum &&
                    TecUtilVarIsEnabled(IntegrationValues.XVarNum))
                {
                    WeightingVar = TecUtilDataValueGetRef(ZoneNum, IntegrationValues.XVarNum);
                }
                else
                {
                    WeightingVar = NULL;
                }
            }
            else
            {
                XComponent = NULL;
                YComponent = NULL;
                ZComponent = NULL;
                WeightingVar = NULL;
            }

            if (UsesScalarVariable)
            {
                ScalarVar = TecUtilDataValueGetRef(ZoneNum, IntegrationValues.ScalarVarNum);
            }
            else
            {
                ScalarVar = NULL;
            }

            IMin = IntegrationValues.IRange.Min;
            JMin = IntegrationValues.JRange.Min;
            KMin = IntegrationValues.KRange.Min;
            IMax = IntegrationValues.IRange.Max;
            JMax = IntegrationValues.JRange.Max;
            KMax = IntegrationValues.KRange.Max;
            if (IMin <= 0) IMin += CurZoneInfo->IMax;
            if (JMin <= 0) JMin += CurZoneInfo->JMax;
            if (KMin <= 0) KMin += CurZoneInfo->KMax;
            if (IMax <= 0) IMax += CurZoneInfo->IMax;
            if (JMax <= 0) JMax += CurZoneInfo->JMax;
            if (KMax <= 0) KMax += CurZoneInfo->KMax;

            /* Alias IntegrateOver to equivalent settings for 2D/1D zones. */
            IntegrateOver = IntegrationValues.IntegrateOver;
            if (CurZoneInfo->ZoneType == ZoneType_Ordered)
            {
                if (CurZoneInfo->KMax == 1)
                {
                    if (CurZoneInfo->JMax == 1)
                    {
                        /* I-ordered zone. */
                        if (IntegrateOver == IntegrateOver_Volume ||
                            IntegrateOver == IntegrateOver_JPlanes ||
                            IntegrateOver == IntegrateOver_KPlanes)
                        {
                            IntegrateOver = IntegrateOver_ILines;
                        }
                        else if (IntegrateOver == IntegrateOver_IPlanes ||
                                 IntegrateOver == IntegrateOver_JLines ||
                                 IntegrateOver == IntegrateOver_KLines)
                        {
                            /* Error */
                            TecUtilDialogErrMsg("Integration of I-ordered zones\n"
                                                "must be along I-Lines.");
                            ReturnStatus = ReturnStatus_Error;
                        }
                    }
                    else if (CurZoneInfo->IMax == 1)
                    {
                        /* J-ordered zone. */
                        if (IntegrateOver == IntegrateOver_Volume ||
                            IntegrateOver == IntegrateOver_IPlanes ||
                            IntegrateOver == IntegrateOver_KPlanes)
                        {
                            IntegrateOver = IntegrateOver_JLines;
                        }
                        else if (IntegrateOver == IntegrateOver_JPlanes ||
                                 IntegrateOver == IntegrateOver_ILines ||
                                 IntegrateOver == IntegrateOver_KLines)
                        {
                            /* Error */
                            TecUtilDialogErrMsg("Integration of J-ordered zones\n"
                                                "must be along J-Lines.");
                            ReturnStatus = ReturnStatus_Error;
                        }
                    }
                    else
                    {
                        /* IJ-ordered zone. */
                        if (IntegrateOver == IntegrateOver_Volume)
                        {
                            IntegrateOver = IntegrateOver_KPlanes;
                        }
                        else if (IntegrateOver == IntegrateOver_IPlanes)
                        {
                            IntegrateOver = IntegrateOver_JLines;
                        }
                        else if (IntegrateOver == IntegrateOver_JPlanes)
                        {
                            IntegrateOver = IntegrateOver_ILines;
                        }
                        else if (IntegrateOver == IntegrateOver_KLines)
                        {
                            /* Error */
                            TecUtilDialogErrMsg("Integration of IJ-ordered zones\n"
                                                "cannot be along K-Lines.");
                            ReturnStatus = ReturnStatus_Error;
                        }
                    }
                }
                else if (CurZoneInfo->JMax == 1)
                {
                    if (CurZoneInfo->IMax == 1)
                    {
                        /* K-ordered zone. */
                        if (IntegrateOver == IntegrateOver_Volume ||
                            IntegrateOver == IntegrateOver_IPlanes ||
                            IntegrateOver == IntegrateOver_JPlanes)
                        {
                            IntegrateOver = IntegrateOver_KLines;
                        }
                        else if (IntegrateOver == IntegrateOver_KPlanes ||
                                 IntegrateOver == IntegrateOver_ILines ||
                                 IntegrateOver == IntegrateOver_JLines)
                        {
                            /* Error */
                            TecUtilDialogErrMsg("Integration of K-ordered zones\n"
                                                "must be along K-Lines.");
                            ReturnStatus = ReturnStatus_Error;
                        }
                    }
                    else
                    {
                        /* IK-ordered zone. */
                        if (IntegrateOver == IntegrateOver_Volume)
                        {
                            IntegrateOver = IntegrateOver_JPlanes;
                        }
                        else if (IntegrateOver == IntegrateOver_IPlanes)
                        {
                            IntegrateOver = IntegrateOver_KLines;
                        }
                        else if (IntegrateOver == IntegrateOver_KPlanes)
                        {
                            IntegrateOver = IntegrateOver_ILines;
                        }
                        else if (IntegrateOver == IntegrateOver_JLines)
                        {
                            /* Error */
                            TecUtilDialogErrMsg("Integration of IK-ordered zones\n"
                                                "cannot be along J-Lines.");
                            ReturnStatus = ReturnStatus_Error;
                        }
                    }
                }
                else if (CurZoneInfo->IMax == 1)
                {
                    /* JK-Ordered zone. */
                    if (IntegrateOver == IntegrateOver_Volume)
                    {
                        IntegrateOver = IntegrateOver_IPlanes;
                    }
                    else if (IntegrateOver == IntegrateOver_JPlanes)
                    {
                        IntegrateOver = IntegrateOver_KLines;
                    }
                    else if (IntegrateOver == IntegrateOver_KPlanes)
                    {
                        IntegrateOver = IntegrateOver_JLines;
                    }
                    else if (IntegrateOver == IntegrateOver_ILines)
                    {
                        /* Error */
                        TecUtilDialogErrMsg("Integration of JK-ordered zones\n"
                                            "cannot be along I-Lines.");
                        ReturnStatus = ReturnStatus_Error;
                    }
                }
            }
            else if (CurZoneInfo->ZoneType == ZoneType_FETriangle ||
                     CurZoneInfo->ZoneType == ZoneType_FEQuad)
            {
                if (IntegrateOver == IntegrateOver_Volume ||
                    IntegrateOver == IntegrateOver_KPlanes)
                {
                    IntegrateOver = IntegrateOver_KPlanes;
                }
                else
                {
                    TecUtilDialogErrMsg("Integration of finite-element zones\n"
                                        "must be over Cell Volumes or K Planes.");
                    ReturnStatus = ReturnStatus_Error;
                }
            }
            else if (CurZoneInfo->ZoneType == ZoneType_FELineSeg)
            {
                if (IntegrateOver == IntegrateOver_Volume ||
                    IntegrateOver == IntegrateOver_JPlanes ||
                    IntegrateOver == IntegrateOver_KPlanes)
                {
                    IntegrateOver = IntegrateOver_ILines;
                }
                else if (IntegrateOver == IntegrateOver_IPlanes ||
                         IntegrateOver == IntegrateOver_JLines ||
                         IntegrateOver == IntegrateOver_KLines)
                {
                    /* Error */
                    TecUtilDialogErrMsg("Integration of FELineSeg zones\n"
                                        "must be along I-Lines.");
                    ReturnStatus = ReturnStatus_Error;
                }
            }
            else
            {
                if (IntegrateOver != IntegrateOver_Volume)
                {
                    TecUtilDialogErrMsg("Integration of 3D finite-element zones\n"
                                        "must be over cell volumes");
                    ReturnStatus = ReturnStatus_Error; /* Only volume integrations allowed for 3D unstructured. */
                }
            }
        }

        /* Disallow axisymmetric volume calculations. */
        if (ReturnStatus == ReturnStatus_OK &&
            IntegrateOver == 0 && IntegrationValues.Axisymmetric)
        {
            TecUtilDialogErrMsg("Cannot perform axisymmetric integrations on 3D zones.");
            ReturnStatus = ReturnStatus_Error;
        }

        if (ReturnStatus == ReturnStatus_OK)
        {
            /* Break the integration assignment for this zone into individual integrations and perform
               them one at a time. */
            sprintf(PercentDialogText, "Integrating zone %d", ZoneNum);
            TecUtilDialogLaunchPercentDone(PercentDialogText, TRUE);
            ResetPercentDone();
            switch (IntegrateOver)
            {
                case IntegrateOver_Volume: /* Volume, can be unstructured */
                    NumberOfIntegrations = 1;
                    ReturnStatus = PerformSingleIntegration(&Numerator,
                                                            &Denominator,
                                                            ScalarVar,
                                                            WeightingVar,
                                                            XComponent,
                                                            YComponent,
                                                            ZComponent,
                                                            IntegrateOver,
                                                            IMin,
                                                            JMin,
                                                            KMin,
                                                            IMax,
                                                            JMax,
                                                            KMax,
                                                            0,
                                                            100);
                    if (ReturnStatus == ReturnStatus_OK)
                    {
                        ZoneNumeratorTotal   = Numerator;
                        ZoneDenominatorTotal = Denominator;

                        if (fabs(Denominator) > FLT_MIN)
                        {
                            Numerator /= Denominator;
                        }

                        if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                                       ZoneNum,
                                                       (LgIndex_t)1,
                                                       (LgIndex_t)1,
                                                       (LgIndex_t)1,
                                                       Numerator))
                        {
                            ReturnStatus = ReturnStatus_Error;
                        }
                    }
                    break;
                case IntegrateOver_IPlanes:
                    NumberOfIntegrations = (IMax - IMin) / IntegrationValues.IRange.Skip + 1;
                    IntegrationNumber = 0;
                    for (I = IMin; ReturnStatus == ReturnStatus_OK && I <= IMax; I += IntegrationValues.IRange.Skip)
                    {
                        IntegrationNumber++;
                        PercentDoneMinimum = 100 * (IntegrationNumber - 1) / NumberOfIntegrations;
                        PercentDoneMaximum = 100 * IntegrationNumber / NumberOfIntegrations;
                        ReturnStatus = PerformSingleIntegration(&Numerator,
                                                                &Denominator,
                                                                ScalarVar,
                                                                WeightingVar,
                                                                XComponent,
                                                                YComponent,
                                                                ZComponent,
                                                                IntegrateOver,
                                                                I,
                                                                JMin,
                                                                KMin,
                                                                I,
                                                                JMax,
                                                                KMax,
                                                                PercentDoneMinimum,
                                                                PercentDoneMaximum);
                        if (ReturnStatus == ReturnStatus_OK)
                        {
                            ZoneNumeratorTotal   += Numerator;
                            ZoneDenominatorTotal += Denominator;

                            if (fabs(Denominator) > FLT_MIN)
                            {
                                Numerator /= Denominator;
                            }

                            if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                                           ZoneNum,
                                                           (LgIndex_t)I,
                                                           (LgIndex_t)1,
                                                           (LgIndex_t)1,
                                                           Numerator))
                            {
                                ReturnStatus = ReturnStatus_Error;
                            }
                        }
                    }
                    break;
                case IntegrateOver_JPlanes:
                    NumberOfIntegrations = (JMax - JMin) / IntegrationValues.JRange.Skip + 1;
                    IntegrationNumber = 0;
                    for (J = JMin; ReturnStatus == ReturnStatus_OK && J <= JMax; J += IntegrationValues.JRange.Skip)
                    {
                        IntegrationNumber++;
                        PercentDoneMinimum = 100 * (IntegrationNumber - 1) / NumberOfIntegrations;
                        PercentDoneMaximum = 100 * IntegrationNumber / NumberOfIntegrations;
                        ReturnStatus = PerformSingleIntegration(&Numerator,
                                                                &Denominator,
                                                                ScalarVar,
                                                                WeightingVar,
                                                                XComponent,
                                                                YComponent,
                                                                ZComponent,
                                                                IntegrateOver,
                                                                IMin,
                                                                J,
                                                                KMin,
                                                                IMax,
                                                                J,
                                                                KMax,
                                                                PercentDoneMinimum,
                                                                PercentDoneMaximum);
                        if (ReturnStatus == ReturnStatus_OK)
                        {
                            ZoneNumeratorTotal   += Numerator;
                            ZoneDenominatorTotal += Denominator;

                            if (fabs(Denominator) > FLT_MIN)
                            {
                                Numerator /= Denominator;
                            }

                            if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                                           ZoneNum,
                                                           (LgIndex_t)1,
                                                           (LgIndex_t)J,
                                                           (LgIndex_t)1,
                                                           Numerator))
                            {
                                ReturnStatus = ReturnStatus_Error;
                            }
                        }
                    }
                    break;
                case IntegrateOver_KPlanes: /* Can be unstructured */
                    IntegrationNumber = 0;
                    if (CurZoneInfo->ZoneType == ZoneType_Ordered)
                    {
                        NumberOfIntegrations = (KMax - KMin) / IntegrationValues.KRange.Skip + 1;
                        for (K = KMin; ReturnStatus == ReturnStatus_OK && K <= KMax; K += IntegrationValues.KRange.Skip)
                        {
                            IntegrationNumber++;
                            PercentDoneMinimum = 100 * (IntegrationNumber - 1) / NumberOfIntegrations;
                            PercentDoneMaximum = 100 * IntegrationNumber / NumberOfIntegrations;
                            ReturnStatus = PerformSingleIntegration(&Numerator,
                                                                    &Denominator,
                                                                    ScalarVar,
                                                                    WeightingVar,
                                                                    XComponent,
                                                                    YComponent,
                                                                    ZComponent,
                                                                    IntegrateOver,
                                                                    IMin,
                                                                    JMin,
                                                                    K,
                                                                    IMax,
                                                                    JMax,
                                                                    K,
                                                                    PercentDoneMinimum,
                                                                    PercentDoneMaximum);
                            if (ReturnStatus == ReturnStatus_OK)
                            {
                                ZoneNumeratorTotal   += Numerator;
                                ZoneDenominatorTotal += Denominator;

                                if (fabs(Denominator) > FLT_MIN)
                                {
                                    Numerator /= Denominator;
                                }

                                if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                                               ZoneNum,
                                                               (LgIndex_t)1,
                                                               (LgIndex_t)1,
                                                               (LgIndex_t)K,
                                                               Numerator))
                                {
                                    ReturnStatus = ReturnStatus_Error;
                                }
                            }
                        }
                    }
                    else
                    {
                        NumberOfIntegrations = 1;
                        ReturnStatus = PerformSingleIntegration(&Numerator,
                                                                &Denominator,
                                                                ScalarVar,
                                                                WeightingVar,
                                                                XComponent,
                                                                YComponent,
                                                                ZComponent,
                                                                IntegrateOver,
                                                                IMin,
                                                                JMin,
                                                                KMin,
                                                                IMax,
                                                                JMax,
                                                                KMin,
                                                                0,
                                                                100);
                        if (ReturnStatus == ReturnStatus_OK)
                        {
                            ZoneNumeratorTotal += Numerator;
                            ZoneDenominatorTotal += Denominator;

                            if (fabs(Denominator) > FLT_MIN)
                            {
                                Numerator /= Denominator;
                            }

                            if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                                           ZoneNum,
                                                           (LgIndex_t)1,
                                                           (LgIndex_t)1,
                                                           (LgIndex_t)1,
                                                           Numerator))
                            {
                                ReturnStatus = ReturnStatus_Error;
                            }
                        }
                    }
                    break;
                case IntegrateOver_ILines:
                    if (CurZoneInfo->ZoneType == ZoneType_Ordered)
                    {
                        NumberOfIntegrations = ((KMax - KMin) / IntegrationValues.KRange.Skip + 1) *
                                               ((JMax - JMin) / IntegrationValues.JRange.Skip + 1);
                        IntegrationNumber = 0;
                        for (J = JMin; ReturnStatus == ReturnStatus_OK && J <= JMax; J += IntegrationValues.JRange.Skip)
                        {
                            for (K = KMin; ReturnStatus == ReturnStatus_OK && K <= KMax; K += IntegrationValues.KRange.Skip)
                            {
                                IntegrationNumber++;
                                PercentDoneMinimum = 100 * (IntegrationNumber - 1) / NumberOfIntegrations;
                                PercentDoneMaximum = 100 * IntegrationNumber / NumberOfIntegrations;
                                ReturnStatus = PerformSingleIntegration(&Numerator,
                                                                        &Denominator,
                                                                        ScalarVar,
                                                                        WeightingVar,
                                                                        XComponent,
                                                                        YComponent,
                                                                        ZComponent,
                                                                        IntegrateOver,
                                                                        IMin,
                                                                        J,
                                                                        K,
                                                                        IMax,
                                                                        J,
                                                                        K,
                                                                        PercentDoneMinimum,
                                                                        PercentDoneMaximum);
                                if (ReturnStatus == ReturnStatus_OK)
                                {
                                    ZoneNumeratorTotal   += Numerator;
                                    ZoneDenominatorTotal += Denominator;

                                    if (fabs(Denominator) > FLT_MIN)
                                    {
                                        Numerator /= Denominator;
                                    }

                                    if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                                                   ZoneNum,
                                                                   (LgIndex_t)1,
                                                                   (LgIndex_t)J,
                                                                   (LgIndex_t)K,
                                                                   Numerator))
                                    {
                                        ReturnStatus = ReturnStatus_Error;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        /* FELineSeg zone */
                        ReturnStatus = PerformSingleIntegration(&Numerator,
                                                                &Denominator,
                                                                ScalarVar,
                                                                WeightingVar,
                                                                XComponent,
                                                                YComponent,
                                                                ZComponent,
                                                                IntegrateOver,
                                                                IMin,
                                                                JMin,
                                                                KMin,
                                                                IMin,
                                                                JMax,
                                                                KMin,
                                                                0,
                                                                100);
                        if (ReturnStatus == ReturnStatus_OK)
                        {
                            ZoneNumeratorTotal   += Numerator;
                            ZoneDenominatorTotal += Denominator;

                            if (fabs(Denominator) > FLT_MIN)
                            {
                                Numerator /= Denominator;
                            }

                            if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                                           ZoneNum,
                                                           (LgIndex_t)1,
                                                           (LgIndex_t)1,
                                                           (LgIndex_t)1,
                                                           Numerator))
                            {
                                ReturnStatus = ReturnStatus_Error;
                            }
                        }
                    }
                    break;
                case IntegrateOver_JLines: /* J lines */
                    NumberOfIntegrations = ((KMax - KMin) / IntegrationValues.KRange.Skip + 1) *
                                           ((IMax - IMin) / IntegrationValues.IRange.Skip + 1);
                    IntegrationNumber = 0;
                    for (I = IMin; ReturnStatus == ReturnStatus_OK && I <= IMax; I += IntegrationValues.IRange.Skip)
                    {
                        for (K = KMin; ReturnStatus == ReturnStatus_OK && K <= KMax; K += IntegrationValues.KRange.Skip)
                        {
                            IntegrationNumber++;
                            PercentDoneMinimum = 100 * (IntegrationNumber - 1) / NumberOfIntegrations;
                            PercentDoneMaximum = 100 * IntegrationNumber / NumberOfIntegrations;
                            ReturnStatus = PerformSingleIntegration(&Numerator,
                                                                    &Denominator,
                                                                    ScalarVar,
                                                                    WeightingVar,
                                                                    XComponent,
                                                                    YComponent,
                                                                    ZComponent,
                                                                    IntegrateOver,
                                                                    I,
                                                                    JMin,
                                                                    K,
                                                                    I,
                                                                    JMax,
                                                                    K,
                                                                    PercentDoneMinimum,
                                                                    PercentDoneMaximum);
                            if (ReturnStatus == ReturnStatus_OK)
                            {
                                ZoneNumeratorTotal   += Numerator;
                                ZoneDenominatorTotal += Denominator;

                                if (fabs(Denominator) > FLT_MIN)
                                {
                                    Numerator /= Denominator;
                                }

                                if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                                               ZoneNum,
                                                               (LgIndex_t)I,
                                                               (LgIndex_t)1,
                                                               (LgIndex_t)K,
                                                               Numerator))
                                {
                                    ReturnStatus = ReturnStatus_Error;
                                }
                            }
                        }
                    }
                    break;
                case IntegrateOver_KLines: /* K lines */
                    NumberOfIntegrations = ((JMax - JMin) / IntegrationValues.JRange.Skip + 1) *
                                           ((IMax - IMin) / IntegrationValues.IRange.Skip + 1);
                    IntegrationNumber = 0;
                    for (I = IMin; I <= IMax; I += IntegrationValues.IRange.Skip)
                    {
                        for (J = JMin; J <= JMax; J += IntegrationValues.JRange.Skip)
                        {
                            IntegrationNumber++;
                            PercentDoneMinimum = 100 * (IntegrationNumber - 1) / NumberOfIntegrations;
                            PercentDoneMaximum = 100 * IntegrationNumber / NumberOfIntegrations;
                            ReturnStatus = PerformSingleIntegration(&Numerator,
                                                                    &Denominator,
                                                                    ScalarVar,
                                                                    WeightingVar,
                                                                    XComponent,
                                                                    YComponent,
                                                                    ZComponent,
                                                                    IntegrateOver,
                                                                    I,
                                                                    J,
                                                                    KMin,
                                                                    I,
                                                                    J,
                                                                    KMax,
                                                                    PercentDoneMinimum,
                                                                    PercentDoneMaximum);
                            if (ReturnStatus == ReturnStatus_OK)
                            {
                                ZoneNumeratorTotal   += Numerator;
                                ZoneDenominatorTotal += Denominator;

                                if (fabs(Denominator) > FLT_MIN)
                                {
                                    Numerator /= Denominator;
                                }

                                if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                                               ZoneNum,
                                                               (LgIndex_t)I,
                                                               (LgIndex_t)J,
                                                               (LgIndex_t)1,
                                                               Numerator))
                                {
                                    ReturnStatus = ReturnStatus_Error;
                                }
                            }
                        }
                    }
                    break;
                default:
                    TecUtilDialogErrMsg("Unrecognized 'Integrate Over' option (internal error).");
                    ReturnStatus = ReturnStatus_Error;
                    break;
            } /* case statement. */
            TecUtilDialogDropPercentDone();
        }

        if (ReturnStatus == ReturnStatus_OK)
        {
            /* Sum all zone totals. */
            NumeratorTotal   += ZoneNumeratorTotal;
            DenominatorTotal += ZoneDenominatorTotal;

            if (NumberOfIntegrations > 1)
            {
                if (fabs(ZoneDenominatorTotal) > FLT_MIN)
                {
                    ZoneNumeratorTotal /= ZoneDenominatorTotal;
                }

                if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                               ZoneNum,
                                               (LgIndex_t)TECINT_TOTAL,
                                               (LgIndex_t)1,
                                               (LgIndex_t)1,
                                               ZoneNumeratorTotal))
                {
                    ReturnStatus = ReturnStatus_Error;
                }
            }
        }
        ZoneNum = (EntIndex_t)TecUtilSetGetNextMember(IntegrationValues.ZoneSet,
                                                      (SetIndex_t)ZoneNum);
    } /* Loop over all zones. */

    if (ReturnStatus == ReturnStatus_OK)
    {
        /* Calculate the weighted average. */
        if (fabs(DenominatorTotal) > FLT_MIN)
        {
            NumeratorTotal /= DenominatorTotal;
        }

        if (!IntegrationResultAddValue(IntegrationValues.IntegrationResult,
                                       (EntIndex_t)TECINT_TOTAL,
                                       (LgIndex_t)1,
                                       (LgIndex_t)1,
                                       (LgIndex_t)1,
                                       NumeratorTotal))
        {
            ReturnStatus = ReturnStatus_Error;
        }
    }

    if (ReturnStatus == ReturnStatus_OK)
    {
        if (ContainsNegativeVolumes && ContainsPositiveVolumes && !TecUtilMacroIsBatchModeActive())
        {
            TecUtilDialogMessageBox("Domain of integration contains both\n"
                                    "positive and negative areas/volumes.",
                                    MessageBox_Warning);
        }
    }

    ENSURE(VALID_ENUM(ReturnStatus, ReturnStatus_e));
    return(ReturnStatus);
}

static Boolean_t ZonesAreValid(Set_pa ZoneSet)
/*
 * If all zones in ZoneSet are valid zone numbers in the current dataset, return TRUE.
 * Otherwise, return FALSE.
 */
{
    SetIndex_t SetCount;
    SetIndex_t SetIndex;
    EntIndex_t ZoneNum;
    Boolean_t  Result = TRUE;

    REQUIRE(VALID_REF(ZoneSet));

    SetCount = TecUtilSetGetMemberCount(ZoneSet);
    for (SetIndex = 1; SetIndex <= SetCount; SetIndex++)
    {
        ZoneNum = (EntIndex_t)TecUtilSetGetMember(ZoneSet, SetIndex);
        if (!TecUtilZoneIsEnabled(ZoneNum))
        {
            Result = FALSE;
            break;
        }
    }

    ENSURE(VALID_BOOLEAN(Result));

    return(Result);
}

int GetZoneDimensionality(EntIndex_t ZoneNum)
{
    ZoneType_e ZoneType;
    int        Dimensionality;

    REQUIRE(TecUtilDataSetIsAvailable());
    REQUIRE(TecUtilZoneIsEnabled(ZoneNum));

    ZoneType = TecUtilZoneGetType(ZoneNum);

    if (ZoneType == ZoneType_Ordered)
    {
        LgIndex_t    IMax;
        LgIndex_t    JMax;
        LgIndex_t    KMax;

        TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax,
                           NULL, NULL, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL);

        Dimensionality = 0;
        if (IMax > 1) Dimensionality++;
        if (JMax > 1) Dimensionality++;
        if (KMax > 1) Dimensionality++;
    }
    else if (ZoneType == ZoneType_FETetra ||
             ZoneType == ZoneType_FEBrick)
    {
        Dimensionality = 3;
    }
    else if (ZoneType == ZoneType_FELineSeg)
    {
        Dimensionality = 1;
    }
    else
    {
        Dimensionality = 2;
    }

    ENSURE(Dimensionality >= 0);
    return(Dimensionality);
}

ReturnStatus_e Integrate(IntegrationResult_pa  IntegrationResult,
                         VariableOption_e      VariableOption,
                         Boolean_t             Axisymmetric,
                         SymmetryVar_e         SymmetryVar,
                         double                SymmetryValue,
                         EntIndex_t            ScalarVarNum,
                         EntIndex_t            XVarNum,
                         EntIndex_t            YVarNum,
                         EntIndex_t            ZVarNum,
                         IntegrateOver_e       IntegrateOver,
                         Set_pa                ZoneSet,
                         IndexRange_t          IRange,
                         IndexRange_t          JRange,
                         IndexRange_t          KRange,
                         Boolean_t             Absolute,
                         Boolean_t             ExcludeBlanked)
{
    char              complaint[1024];
    char              EnvironmentVariable[1024];
    double            Result;
    EntIndex_t        NumZones;
    EntIndex_t        NumVars;
    EntIndex_t        ZoneNum;
    LgIndex_t         IMin;
    LgIndex_t         JMin;
    LgIndex_t         KMin;
    LgIndex_t         IMax;
    LgIndex_t         JMax;
    LgIndex_t         KMax;
    LgIndex_t         ZoneIMax;
    LgIndex_t         ZoneJMax;
    LgIndex_t         ZoneKMax;
    ZoneType_e        ZoneType;
    FrameMode_e       FrameMode = TecUtilFrameGetMode();
    Boolean_t         IsVectorIntegral;
    Boolean_t         AllocatedZoneSet = FALSE;
    ReturnStatus_e    ReturnStatus = ReturnStatus_OK;

    IntegrationValues.IntegrationResult = IntegrationResult;
    IntegrationValues.VariableOption    = VariableOption;
    IntegrationValues.Axisymmetric      = Axisymmetric;
    IntegrationValues.SymmetryVar       = SymmetryVar;
    IntegrationValues.SymmetryValue     = SymmetryValue;
    IntegrationValues.ScalarVarNum      = ScalarVarNum;
    IntegrationValues.XVarNum           = XVarNum;
    IntegrationValues.YVarNum           = YVarNum;
    IntegrationValues.ZVarNum           = ZVarNum;
    IntegrationValues.IntegrateOver     = IntegrateOver;
    IntegrationValues.ZoneSet           = ZoneSet;
    IntegrationValues.IRange            = IRange;
    IntegrationValues.JRange            = JRange;
    IntegrationValues.KRange            = KRange;
    IntegrationValues.Absolute          = Absolute;
    IntegrationValues.ExcludeBlanked    = ExcludeBlanked;

    IsVectorIntegral = (IntegrationValues.VariableOption == VariableOption_VectorDotNormal ||
                        IntegrationValues.VariableOption == VariableOption_VectorAverage ||
                        IntegrationValues.VariableOption == VariableOption_VectorDotTangential);

    if (ReturnStatus == ReturnStatus_OK)
    {
        /* Check zone list for integration. */
        if (!TecUtilDataSetIsAvailable())
        {
            TecUtilDialogErrMsg("No data set available in the current frame.");
            ReturnStatus = ReturnStatus_Error;
        }
        else
        {
            TecUtilDataSetGetInfo(NULL,
                                  &NumZones,
                                  &NumVars);
        }
    }

    /* Check input parameters for obvious errors likely due to bad programming (as opposed to user input). */
    if (!VALID_REF(IntegrationValues.IntegrationResult))
    {
        ReturnStatus = ReturnStatus_BadIntegrationResult;
    }
    else if (!VALID_ENUM(IntegrationValues.VariableOption, VariableOption_e))
    {
        ReturnStatus = ReturnStatus_BadVariableOption;
    }
    else if (!VALID_BOOLEAN(IntegrationValues.Axisymmetric))
    {
        ReturnStatus = ReturnStatus_BadAxisymmetric;
    }
    else if (IntegrationValues.Axisymmetric &&
             !VALID_ENUM(IntegrationValues.SymmetryVar, SymmetryVar_e))
    {
        ReturnStatus = ReturnStatus_BadSymmetryVar;
    }
    else if ((IntegrationValues.VariableOption == VariableOption_ScalarIntegral ||
              IntegrationValues.VariableOption == VariableOption_Average ||
              IntegrationValues.VariableOption == VariableOption_WeightedAverage) &&
             (IntegrationValues.ScalarVarNum <= 0 ||
              !TecUtilVarIsEnabled(IntegrationValues.ScalarVarNum)))
    {
        ReturnStatus = ReturnStatus_BadScalarVarNum;
    }
    else if ((IntegrationValues.VariableOption == VariableOption_WeightedAverage ||
              IsVectorIntegral) &&
             (IntegrationValues.XVarNum <= 0 ||
              !TecUtilVarIsEnabled(IntegrationValues.XVarNum)))
    {
        ReturnStatus = ReturnStatus_BadXVarNum;
    }
    else if (IsVectorIntegral &&
             (IntegrationValues.YVarNum <= 0 ||
              !TecUtilVarIsEnabled(IntegrationValues.YVarNum)))
    {
        ReturnStatus = ReturnStatus_BadYVarNum;
    }
    else if (IsVectorIntegral &&
             FrameMode == Frame_ThreeD &&
             (IntegrationValues.ZVarNum <= 0 ||
              !TecUtilVarIsEnabled(IntegrationValues.ZVarNum)))
    {
        ReturnStatus = ReturnStatus_BadZVarNum;
    }
    else if (!VALID_ENUM(IntegrationValues.IntegrateOver, IntegrateOver_e))
    {
        ReturnStatus = ReturnStatus_BadIntegrateOver;
    }
    else if (!VALID_REF(IntegrationValues.ZoneSet) && IntegrationValues.ZoneSet != NULL)
    {
        ReturnStatus = ReturnStatus_BadZoneSet;
    }
    else if (!VALID_BOOLEAN(IntegrationValues.Absolute))
    {
        ReturnStatus = ReturnStatus_BadAbsolute;
    }
    else if (!VALID_BOOLEAN(IntegrationValues.ExcludeBlanked))
    {
        ReturnStatus = ReturnStatus_BadExcludeBlanked;
    }

    if (ReturnStatus == ReturnStatus_OK)
    {
        /* Make sure we are in an appropriate frame mode. */
        if (FrameMode == Frame_Sketch)
        {
            TecUtilDialogErrMsg("Cannot integrate while in Sketch frame mode.");
            ReturnStatus = ReturnStatus_Error;
        }
        else if (FrameMode == Frame_XY)
        {
            /* Cannot do vector calculations in XY frame mode. */
            if (IsVectorIntegral)
            {
                TecUtilDialogErrMsg("Cannot perform vector-type integrations in XY frame mode.");
                ReturnStatus = ReturnStatus_Error;
            }
        }
    }

    if (ReturnStatus == ReturnStatus_OK && !IntegrationValues.ZoneSet)
    {
        /* Integrate over all valid (enabled) zones. */
        if (TecUtilZoneGetEnabled(&IntegrationValues.ZoneSet))
        {
            AllocatedZoneSet = TRUE;
        }
        else
        {
            TecUtilDialogErrMsg("Error getting enabled zone set.");
            ReturnStatus = ReturnStatus_Error;
        }
    }

    if (ReturnStatus == ReturnStatus_OK && !ZonesAreValid(IntegrationValues.ZoneSet))
    {
        TecUtilDialogErrMsg("The integration zone set contains invalid zone numbers.");
        ReturnStatus = ReturnStatus_BadZoneSet;
    }

    if (ReturnStatus == ReturnStatus_OK)
    {
        /* Make sure Axisymmetric isn't set for a 3D frame mode. */
        if (IntegrationValues.Axisymmetric &&
            TecUtilFrameGetMode() == Frame_ThreeD)
        {
            TecUtilDialogErrMsg("Axisymmetric integration cannot be performed\n"
                                "when the frame mode is 3D.");
            ReturnStatus = ReturnStatus_Error;
        }
    }

    if (ReturnStatus == ReturnStatus_OK)
    {
        /* Check skip values. */
        if (IntegrationValues.IRange.Skip < 1)
        {
            TecUtilDialogErrMsg("Invalid I Skip value");
            ReturnStatus = ReturnStatus_BadIRange;
        }
        else if (IntegrationValues.JRange.Skip < 1)
        {
            TecUtilDialogErrMsg("Invalid J Skip value");
            ReturnStatus = ReturnStatus_BadJRange;
        }
        else if (IntegrationValues.KRange.Skip < 1)
        {
            TecUtilDialogErrMsg("Invalid K Skip value");
            ReturnStatus = ReturnStatus_BadKRange;
        }
    }

    if (ReturnStatus == ReturnStatus_OK)
    {
        /* Check index limits and integration type for all zones. */
        ZoneNum = (EntIndex_t)TecUtilSetGetNextMember(IntegrationValues.ZoneSet,
                                                      TECUTILSETNOTMEMBER);
        while (ReturnStatus == ReturnStatus_OK && ZoneNum != (EntIndex_t)TECUTILSETNOTMEMBER)
        {
            int Dimensionality;

            if (!TecUtilZoneIsEnabled(ZoneNum))
            {
                sprintf(complaint, "Zone %d is not enabled. Cannot perform integration.", ZoneNum);
                TecUtilDialogErrMsg(complaint);
                ReturnStatus = ReturnStatus_Error;
            }

            if (ReturnStatus == ReturnStatus_OK)
            {
                TecUtilZoneGetInfo(ZoneNum, &ZoneIMax, &ZoneJMax, &ZoneKMax,
                                   NULL, NULL, NULL, NULL, NULL,
                                   NULL, NULL, NULL, NULL, NULL);
                Dimensionality = GetZoneDimensionality(ZoneNum);
                ZoneType = TecUtilZoneGetType(ZoneNum);
                JMin = (IntegrationValues.JRange.Min > 0 ? IntegrationValues.JRange.Min : ZoneJMax + IntegrationValues.JRange.Min);
                JMax = (IntegrationValues.JRange.Max > 0 ? IntegrationValues.JRange.Max : ZoneJMax + IntegrationValues.JRange.Max);
                if (ZoneType == ZoneType_Ordered)
                {
                    IMin = (IntegrationValues.IRange.Min > 0 ? IntegrationValues.IRange.Min : ZoneIMax + IntegrationValues.IRange.Min);
                    KMin = (IntegrationValues.KRange.Min > 0 ? IntegrationValues.KRange.Min : ZoneKMax + IntegrationValues.KRange.Min);
                    IMax = (IntegrationValues.IRange.Max > 0 ? IntegrationValues.IRange.Max : ZoneIMax + IntegrationValues.IRange.Max);
                    KMax = (IntegrationValues.KRange.Max > 0 ? IntegrationValues.KRange.Max : ZoneKMax + IntegrationValues.KRange.Max);
                }
                else
                {
                    IMin = 1;
                    KMin = 1;
                    IMax = 1;
                    KMax = 1;
                    ZoneIMax = 1;
                    ZoneKMax = 1;
                }

                if (IMin < 1 ||
                    JMin < 1 ||
                    KMin < 1 ||
                    IMax < 1 ||
                    JMax < 1 ||
                    KMax < 1)
                {
                    sprintf(complaint, "Negative index limits detected for zone %d", ZoneNum);
                    TecUtilDialogErrMsg(complaint);
                    ReturnStatus = ReturnStatus_Error;
                }
            }

            if (ReturnStatus == ReturnStatus_OK && IMax > ZoneIMax)
            {
                sprintf(complaint, "IMax > Zone IMax for zone %d", ZoneNum);
                TecUtilDialogErrMsg(complaint);
                ReturnStatus = ReturnStatus_Error;
            }

            if (ReturnStatus == ReturnStatus_OK && JMax > ZoneJMax)
            {
                sprintf(complaint, "JMax > Zone JMax for zone %d", ZoneNum);
                TecUtilDialogErrMsg(complaint);
                ReturnStatus = ReturnStatus_Error;
            }

            if (ReturnStatus == ReturnStatus_OK && KMax > ZoneKMax)
            {
                sprintf(complaint, "KMax > Zone KMax for zone %d", ZoneNum);
                TecUtilDialogErrMsg(complaint);
                ReturnStatus = ReturnStatus_Error;
            }

            if (ReturnStatus == ReturnStatus_OK && IMin > IMax)
            {
                sprintf(complaint, "IMin > IMax for zone %d", ZoneNum);
                TecUtilDialogErrMsg(complaint);
                ReturnStatus = ReturnStatus_Error;
            }

            if (ReturnStatus == ReturnStatus_OK && JMin > JMax)
            {
                sprintf(complaint, "JMin > JMax for zone %d", ZoneNum);
                TecUtilDialogErrMsg(complaint);
                ReturnStatus = ReturnStatus_Error;
            }

            if (ReturnStatus == ReturnStatus_OK && KMin > KMax)
            {
                sprintf(complaint, "KMin > KMax for zone %d", ZoneNum);
                TecUtilDialogErrMsg(complaint);
                ReturnStatus = ReturnStatus_Error;
            }

            /* I/J/K planes have no tangential direction in 3D (in 2D they become lines, which is ok). */
            if (ReturnStatus == ReturnStatus_OK &&
                ZoneType == ZoneType_Ordered &&
                Dimensionality == 3 &&
                IntegrationValues.VariableOption == VariableOption_VectorDotTangential &&
                (IntegrationValues.IntegrateOver == IntegrateOver_IPlanes ||
                 IntegrationValues.IntegrateOver == IntegrateOver_JPlanes ||
                 IntegrationValues.IntegrateOver == IntegrateOver_KPlanes))
            {
                TecUtilDialogErrMsg("Cannot integrate vector-dot-tangential\n"
                                    "over planes in three dimensions (there\n"
                                    "is no identifiable tangential direction)");
                ReturnStatus = ReturnStatus_Error;
            }

            /* I/J/K planes and lines do not exist for unstructured zones. */
            if (ReturnStatus == ReturnStatus_OK &&
                (ZoneType == ZoneType_FEBrick ||
                 ZoneType == ZoneType_FETetra))
            {
                if (IntegrationValues.IntegrateOver != IntegrateOver_Volume)
                {
                    TecUtilDialogErrMsg("Can only do volume integration for 3D finite-element zones.");
                    ReturnStatus = ReturnStatus_Error;
                }
            }

            if (ReturnStatus == ReturnStatus_OK &&
                (ZoneType == ZoneType_FEQuad ||
                 ZoneType == ZoneType_FETriangle))
            {
                if (IntegrationValues.IntegrateOver != IntegrateOver_Volume &&
                    IntegrationValues.IntegrateOver != IntegrateOver_KPlanes)
                {
                    TecUtilDialogErrMsg("Can only do K-face/volume integration for 2D finite-element zones.");
                    ReturnStatus = ReturnStatus_Error;
                }
            }

            if (ReturnStatus == ReturnStatus_OK &&
                ZoneType == ZoneType_FELineSeg)
            {
                if (IntegrationValues.IntegrateOver == IntegrateOver_IPlanes ||
                    IntegrationValues.IntegrateOver == IntegrateOver_JLines ||
                    IntegrationValues.IntegrateOver == IntegrateOver_KLines)
                {
                    TecUtilDialogErrMsg("Cannot do I-plane, J-line or K-line integration for line segment zones.");
                    ReturnStatus = ReturnStatus_Error;
                }
            }

            ZoneNum = (EntIndex_t)TecUtilSetGetNextMember(IntegrationValues.ZoneSet,
                                                          (SetIndex_t)ZoneNum);
        }
    }

    if (ReturnStatus == ReturnStatus_OK)
    {
        IntegrationResultInitialize(IntegrationValues.IntegrationResult,
                                    IntegrationValues.IntegrateOver);
        switch (DoIntegration())
        {
            case ReturnStatus_Error:
                TecUtilDialogErrMsg("Error occurred during integration.");
                ReturnStatus = ReturnStatus_Error;
                break;
            case ReturnStatus_Canceled:
                TecUtilDialogErrMsg("Integration cancelled.");
                ReturnStatus = ReturnStatus_Error;
                break;
            default:
                break;
        }
    }

    /* Save results to environment variables (maybe eventually to macro variables */
#ifdef WIN32
#define putenv _putenv
#endif
    if (ReturnStatus == ReturnStatus_OK)
    {
        int       ResultCount = IntegrationResultGetCount(IntegrationValues.IntegrationResult);
        LgIndex_t I;
        LgIndex_t J;
        LgIndex_t K;

        if (IntegrationResultGetNumberedValue(IntegrationValues.IntegrationResult,
                                              ResultCount,
                                              &ZoneNum,
                                              &I,
                                              &J,
                                              &K,
                                              &Result))
        {
            sprintf(EnvironmentVariable, "INTEGRATION_TOTAL=%.15g", Result);
            putenv(EnvironmentVariable);
        }
    }

    if (AllocatedZoneSet)
    {
        TecUtilSetDealloc(&IntegrationValues.ZoneSet);
    }

    ENSURE(VALID_ENUM(ReturnStatus, ReturnStatus_e));
    return(ReturnStatus);
}

