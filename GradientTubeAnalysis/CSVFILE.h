#ifndef CSVFILE_H_
#define CSVFILE_H_




Boolean_t CSVCreateFileAndInsertData(char* PathName, std::vector<std::vector<std::vector<ImportType_t> > > *Data,
	std::vector<EntIndex_t> *IntegrationVarNums);

#endif