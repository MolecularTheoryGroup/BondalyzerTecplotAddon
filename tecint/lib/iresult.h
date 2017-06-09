
#if !defined __IRESULT_H
#define __IRESULT_H

extern void IntegrationResultInitialize(IntegrationResult_pa IntegrationResult,
                                        IntegrateOver_e      IntegrateOver);

extern Boolean_t IntegrationResultAddValue(IntegrationResult_pa IntegrationResult,
                                           EntIndex_t           Zone,
                                           LgIndex_t            I,
                                           LgIndex_t            J,
                                           LgIndex_t            K,
                                           double               Value);

#endif /* __IRESULT_H */
