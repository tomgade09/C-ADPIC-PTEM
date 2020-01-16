#include "EField/EField.h"
#include "EField/QSPS.h"
//#include "EField/AlfvenLUT.h"

#include <sstream>
#include <filesystem>

//CUDA includes
#include "device_launch_parameters.h"
#include "ErrorHandling/cudaErrorCheck.h"
#include "ErrorHandling/cudaDeviceMacros.h"

#include "utils/serializationHelpers.h"

using std::to_string;
using std::stringstream;
using namespace utils::fileIO::serialize;

__host__ __device__ EElem::EElem(const char* modelName) : modelName_m{ modelName }
{

}

__host__ __device__ virtual EElem::~EElem()
{

}

__host__ virtual string EElem::name() const
{
	return modelName_m;
}

__host__ virtual EElem** EElem::getPtrGPU() const
{
	return this_d;
}


//EField ctor, dtor
__host__ __device__ EField::EField(bool useGPU) : useGPU_m{ useGPU }
{
	#ifndef __CUDA_ARCH__ //host code
	if(useGPU_m) setupEnvironment();
	#endif /* !__CUDA_ARCH__ */
}

__host__ EField::EField(string serialFolder)
{
	deserialize(serialFolder);
	if(useGPU_m) setupEnvironment();
}

__host__ __device__ EField::~EField()
{
	#ifndef __CUDA_ARCH__ //host code
	if(useGPU_m)deleteEnvironment();
	#endif /* !__CUDA_ARCH__ */
}


//device global kernels
__global__ void setupEnvironmentGPU_EField(EField** efield, EElem*** eelems)
{
	ZEROTH_THREAD_ONLY("setupEnvironmentGPU_EField",
		(*efield) = new EField();
		(*efield)->elemArray(eelems);
	);
}

__global__ void deleteEnvironmentGPU_EField(EField** efield)
{
	ZEROTH_THREAD_ONLY("deleteEnvironmentGPU_EField", delete (*efield));
}

__global__ void addGPU_EField(EField** efield, EElem** elem)
{
	ZEROTH_THREAD_ONLY("addGPU_EField", (*efield)->add(elem));
}

__global__ void increaseCapacity_EField(EField** efield, EElem*** newArray, int capacity)
{
	ZEROTH_THREAD_ONLY("increaseCapacity_EField",
		EElem*** oldArray{ (*efield)->elemArray() };

		for (int elem = 0; elem < (*efield)->size(); elem++)
			newArray[elem] = oldArray[elem];

		(*efield)->capacity(capacity);
		(*efield)->elemArray(newArray); //still retaining the pointer to this memory on host, so no big deal if it's lost here
	);
}


//EField functions
__host__ string EField::getEElemsStr() const
{
	stringstream out;
	for (int elem = 0; elem < size_d; elem++) { out << element(elem)->name() << ", "; }
	return out.str();
}

void EField::setupEnvironment()
{
	CUDA_API_ERRCHK(cudaMalloc((void **)&this_d, sizeof(EField**)));              //allocate memory for EField**
	CUDA_API_ERRCHK(cudaMalloc((void**)&Eelems_d, sizeof(EElem**) * capacity_d)); //allocate memory for EElem** array
	CUDA_API_ERRCHK(cudaMemset(Eelems_d, 0, sizeof(EElem**) * capacity_d));       //clear memory

	setupEnvironmentGPU_EField <<< 1, 1 >>> (this_d, Eelems_d);
	CUDA_KERNEL_ERRCHK_WSYNC();
}

void EField::deleteEnvironment()
{
	deleteEnvironmentGPU_EField <<< 1, 1 >>> (this_d);
	CUDA_KERNEL_ERRCHK_WSYNC();

	CUDA_API_ERRCHK(cudaFree(this_d));
	CUDA_API_ERRCHK(cudaFree(Eelems_d));
}

//#ifndef __CUDA_ARCH__ //host code
__host__ EElem* EField::element(int ind) const
{
	return Eelems_m.at(ind).get();
}

__host__ void EField::add(unique_ptr<EElem> eelem)
{
	if (capacity_d == size_d)
	{
		EElem*** oldArray{ Eelems_d }; //retain so we can cudaFree at the end
		capacity_d += 5;
		
		CUDA_API_ERRCHK(cudaMalloc((void**)&Eelems_d, sizeof(EElem**) * capacity_d)); //create new array that is 5 larger in capacity than the previous
		CUDA_API_ERRCHK(cudaMemset(Eelems_d, 0, sizeof(EElem**) * capacity_d));

		increaseCapacity_EField <<< 1, 1 >>> (this_d, Eelems_d, capacity_d);
		CUDA_KERNEL_ERRCHK();

		CUDA_API_ERRCHK(cudaFree(oldArray));
	}

	//add elem to dev
	addGPU_EField <<< 1, 1 >>> (this_d, eelem->getPtrGPU());
	CUDA_KERNEL_ERRCHK_WSYNC();
	
	//add elem to host
	Eelems_m.push_back(std::move(eelem));
	size_d++;
}
//#endif /* !__CUDA_ARCH__ */

__device__ void EField::add(EElem** newElem)
{
	Eelems_d[size_d] = newElem;
	size_d++;
}

__host__ __device__ int EField::capacity() const
{
	return capacity_d;
}

__host__ __device__ int EField::size() const
{
	return size_d;
}

__device__ void EField::capacity(int cap)
{
	capacity_d = cap;
}

__device__ EElem*** EField::elemArray() const
{
	return Eelems_d;
}

__device__ void EField::elemArray(EElem*** eelems)
{
	Eelems_d = eelems;
}
	
__host__ EField** EField::getPtrGPU() const
{
	return this_d;
}

__host__ __device__ Vperm EField::getEFieldAtS(const meters s, const seconds t) const
{
	tesla ret{ 0.0 };

	#ifndef __CUDA_ARCH__ //host code
	for (auto& elem : Eelems_m) //vector of unique_ptr<EElem>'s
		ret += elem->getEFieldAtS(s, t);
	#else //device code
	for (int elem = 0; elem < size_d; elem++) //c-style array of EElem*'s
		ret += (*(Eelems_d[elem]))->getEFieldAtS(s, t);
	#endif /* !__CUDA_ARCH__ */

	return ret;
}

__host__ void EField::serialize(string serialFolder) const
{
	string filename{ serialFolder + string("/EField.ser") };

	if (std::filesystem::exists(filename))
		cerr << "EField::serialize: Warning: filename exists: " << filename << " You are overwriting an existing file.\n";

	ofstream out(filename, std::ofstream::binary);
	if (!out) throw invalid_argument("EField::serialize: unable to create file: " + filename);
	
	auto writeStrBuf = [&](const stringbuf& sb)
	{
		out.write(sb.str().c_str(), sb.str().length());
	};
	
	size_t size{ Eelems_m };
	out.write(reinterpret_cast<char*>(&size), sizeof(size_t));

	int QSPScnt{ 0 };
	int ALUTcnt{ 0 };
	for (const auto& elem : Eelems_m)
	{
		if (elem.name() == "QSPS") ++QSPScnt;//writeStrBuf(serializeString(string(elem.name()) + to_string(++QSPScnt)));
		else if (elem.name() == "AlfvenLUT") ++ALUTcnt;//writeStrBuf(serializeString(string(elem.name()) + to_string(++ALUTcnt)));
		else throw invalid_argument("EField::serialize: element does not have a recognized name: " + string(elem.name()));
	}

	out.write(reinterpret_cast<char*>(&QSPScnt), sizeof(int));
	out.write(reinterpret_cast<char*>(&ALUTcnt), sizeof(int));

	for (const auto elem : Eelems_m)
	{
		elem.serialize(serialFolder);
		/*  //need to verify use of filename::exists and filename::move below
		if (elem.name() == "QSPS" && std::filename::exists("EField_QSPS.ser"))
		{	int iter{ 0 };
			while (std::filename::exists("EField_QSPS" + to_string(iter) + ".ser")) iter++;
			std::filename::move("EField_QSPS.ser", "EField_QSPS" + to_string(iter) + ".ser");
		}
		else if (elem.name() == "AlfvenLUT" && std::filename::exists("EField_AlfvenLUT.ser"))
		{
			int iter{ 0 };
			while (std::filename::exists("EField_AlfvenLUT" + to_string(iter) + ".ser")) iter++;
			std::filename::move("EField_AlfvenLUT.ser", "EField_AlfvenLUT" + to_string(iter) + ".ser");
		}
		*/
	}

	out.close();
}

__host__ void EField::deserialize(string serialFolder) const
{
	string filename{ serialFolder + string("/EField.ser") };
	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument("EField::deserialize: unable to open file: " + filename);

	size_t elemcnt{ readSizetLength(in) };

	int QSPScnt{ 0 };
	int ALUTcnt{ 0 };

	in.read(&QSPScnt, sizeof(int));
	in.read(&ALUTcnt, sizeof(int));

	for (size_t elem = 0; elem < QSPScnt; elem++)
		Eelems_m.push_back(std::move(std::make_unique<QSPS>(serialFolder, elem)));

	//for (size_t elem = 0; elem < ALUTcnt; elem++)
		//Eelems_m.push_back(std::move(std::make_unique<AlfvenLUT>(serialFolder, elem)));
}