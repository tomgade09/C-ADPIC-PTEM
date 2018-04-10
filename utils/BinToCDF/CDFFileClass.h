#ifndef CDFFILE_H
#define CDFFILE_H

#include <iostream>
#include <string>
#include <vector>

#define WIN32 //for cdf
#include "cdf.h"
#undef WIN32

#define CDFCHKERR() { if (exitStatus_m != CDF_OK) { CDFgetStatusText(exitStatus_m, errTxt_m); std::cout << "Something went wrong: " << exitStatus_m << ":" << errTxt_m << " : " << __FILE__ << ":" << __LINE__ << std::endl; } }

class CDFFileClass //using "Class" just in case there would be any namespace conflicts with the CDF source
{
protected:
	CDFid       cdfid_m;
	std::string filename_m;

	//Attribute and variable characteristics
	std::vector<long>        attrIndicies_m;
	std::vector<std::string> attrNameStrs_m;
	std::vector<long>        zVarIndicies_m;
	std::vector<std::string> zVarNameStrs_m;

	//Error handling
	CDFstatus   exitStatus_m;
	char        errTxt_m[CDF_STATUSTEXT_LEN + 1];

	long findzVarCDFIndexByName(std::string varName);

public:
	CDFFileClass(std::string filename) : filename_m{ filename }
	{
		char fn[1024];
		for (int chr = 0; chr < filename_m.size() + 1; chr++)
			fn[chr] = filename_m.c_str()[chr];
		
		exitStatus_m = CDFcreateCDF(fn, &cdfid_m); //file exists is -2013

		while (exitStatus_m == -2013)
		{
			if (filename_m.size() > 1023)
				throw std::invalid_argument("CDFFileClass::CDFFileClass: filename is too long " + filename);

			filename_m += "2";
			for (int chr = 0; chr < filename_m.size() + 1; chr++)
				fn[chr] = filename_m.c_str()[chr];
			exitStatus_m = CDFcreateCDF(fn, &cdfid_m);
			
		}
	}

	~CDFFileClass()
	{
		CDFcloseCDF(cdfid_m); CDFCHKERR();
	}

	void writeNewZVar(std::string varName, long cdftype, std::vector<int> dimSizes, void* arrayXD);
	void createGlobalAttr(std::string attrName, long cdftype, long datalen, void* data);
	void createZVarAttr(long varNum,         std::string attrName, long cdftype, long datalen, void* data);
	void createZVarAttr(std::string varName, std::string attrName, long cdftype, long datalen, void* data);
};

#endif /* !CDFFILE_H */