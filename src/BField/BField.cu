#include "BField/BField.h"

__host__ __device__ BField::BField(const char* modelName) : name_m{ modelName }
{
    
}

__host__ __device__ virtual BField::~BField()
{

}

__host__ virtual std::string BField::name() const
{
    return modelName_m;
}

__host__ virtual BField** BField::getPtrGPU() const
{
    return this_d;
}