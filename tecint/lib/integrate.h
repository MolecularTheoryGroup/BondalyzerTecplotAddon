static struct IntegrationValues_s
{
    char *           FrameName;
    VariableOption_e VariableOption;
    double           XOrigin;
    double           YOrigin;
    double           ZOrigin;
    EntIndex_t       ScalarVarNum;
    EntIndex_t       XVarNum;
    EntIndex_t       YVarNum;
    EntIndex_t       ZVarNum;
    IntegrateOver_e  IntegrateOver;
    Set_pa           ZoneSet;
    IndexRange_t     IRange;
    IndexRange_t     JRange;
    IndexRange_t     KRange;
    Boolean_t        Absolute;
    Boolean_t        ExcludeBlanked;
    Boolean_t        ShowResults;
    Boolean_t        PlotResults;
    char             PlotAs[MAX_STRING_LENGTH + 1];
    char *           ResultsText;
    char             SaveFile[MAX_PATH_LENGTH + 1];
} IntegrationValues_t;
