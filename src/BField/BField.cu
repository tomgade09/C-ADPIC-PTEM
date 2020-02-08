#include "BField/BField.h"

__host__ __device__ BField::BField(const char* modelName) : name_m{ modelName }
{
    
}

__host__ __device__ BField::~BField()
{

}

__host__ string BField::name() const
{
    return name_m;
}

__host__ BField** BField::getPtrGPU() const
{
    return this_d;
}