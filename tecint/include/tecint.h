
#if !defined __TECINT_H
#define __TECINT_H

#define TECINT_TOTAL (-1)

#ifdef __cplusplus
extern "C"
{
#endif

    typedef struct IntegrationResult_s *IntegrationResult_pa;

    typedef enum
    {
        VariableOption_LengthAreaVolume = 1,
        VariableOption_ScalarIntegral,
        VariableOption_Average,
        VariableOption_WeightedAverage,
        VariableOption_VectorDotNormal,
        VariableOption_VectorAverage,
        VariableOption_VectorDotTangential,
        END_VariableOption_e
    } VariableOption_e;

    typedef enum {
        SymmetryVar_X = 1,
        SymmetryVar_Y,
        END_SymmetryVar_e
    } SymmetryVar_e;

    typedef struct
    {
        int Min;
        int Max;
        int Skip;
    } IndexRange_t;

    typedef enum {
        IntegrateOver_Volume = 1,
        IntegrateOver_IPlanes,
        IntegrateOver_JPlanes,
        IntegrateOver_KPlanes,
        IntegrateOver_ILines,
        IntegrateOver_JLines,
        IntegrateOver_KLines,
        END_IntegrateOver_e
    } IntegrateOver_e;

    typedef enum
    {
        ReturnStatus_OK,
        ReturnStatus_BadIntegrationResult,
        ReturnStatus_BadVariableOption,
        ReturnStatus_BadAxisymmetric,
        ReturnStatus_BadSymmetryVar,
        ReturnStatus_BadSymmetryValue,
        ReturnStatus_BadScalarVarNum,
        ReturnStatus_BadXVarNum,
        ReturnStatus_BadYVarNum,
        ReturnStatus_BadZVarNum,
        ReturnStatus_BadIntegrateOver,
        ReturnStatus_BadZoneSet,
        ReturnStatus_BadIRange,
        ReturnStatus_BadJRange,
        ReturnStatus_BadKRange,
        ReturnStatus_BadAbsolute,
        ReturnStatus_BadExcludeBlanked,
        ReturnStatus_Canceled,
        ReturnStatus_Error,
        END_ReturnStatus_e
    } ReturnStatus_e;

    extern Boolean_t IntegrationResultAlloc(IntegrationResult_pa *IntegrationResult);

    extern void IntegrationResultDealloc(IntegrationResult_pa *IntegrationResult);

    extern int IntegrationResultGetCount(IntegrationResult_pa IntegrationResult);

    extern Boolean_t IntegrationResultGetValue(IntegrationResult_pa  IntegrationResult,
                                               EntIndex_t            Zone,
                                               LgIndex_t             I,
                                               LgIndex_t             J,
                                               LgIndex_t             K,
                                               double               *Value);

    extern Boolean_t IntegrationResultGetFirstValue(IntegrationResult_pa  IntegrationResult,
                                                    EntIndex_t           *Zone,
                                                    LgIndex_t            *I,
                                                    LgIndex_t            *J,
                                                    LgIndex_t            *K,
                                                    double               *Value);

    extern Boolean_t IntegrationResultGetNextValue(IntegrationResult_pa  IntegrationResult,
                                                   EntIndex_t           *Zone,
                                                   LgIndex_t            *I,
                                                   LgIndex_t            *J,
                                                   LgIndex_t            *K,
                                                   double               *Value);

    extern Boolean_t IntegrationResultGetNumberedValue(IntegrationResult_pa  IntegrationResult,
                                                       int                   WhichValue,
                                                       EntIndex_t           *Zone,
                                                       LgIndex_t            *I,
                                                       LgIndex_t            *J,
                                                       LgIndex_t            *K,
                                                       double               *Value);

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
                             Boolean_t             ExcludeBlanked);
#ifdef __cplusplus
}
#endif

#endif /* __TECINT_H */

