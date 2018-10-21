#ifndef __AF_UTILITIES_ENUMERATE_CUH__
#define __AF_UTILITIES_ENUMERATE_CUH__

#include <cuda_runtime.h>

#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include "../particle.cuh"

namespace af
{

    typedef thrust::counting_iterator<uint> counter;
    typedef thrust::tuple<uint, Particle> enumeratedParticle;
    typedef ParticleDeviceArray::iterator deviceParticlesIter;
    typedef thrust::tuple<counter, deviceParticlesIter> deviceEnumeratedParticleIter;
    typedef thrust::zip_iterator<deviceEnumeratedParticleIter> deviceEnumeratedParticlesZipIter;

    class EnumeratedDeviceParticles
    {
        deviceEnumeratedParticlesZipIter _begin;
        deviceEnumeratedParticlesZipIter _end;
    public:
        EnumeratedDeviceParticles(ParticleDeviceArray& particles)
        {
            counter idx_begin(0);
            counter idx_end = idx_begin + particles.size();
            _begin = thrust::make_zip_iterator(thrust::make_tuple(idx_begin, particles.begin()));
            _end = thrust::make_zip_iterator(thrust::make_tuple(idx_end, particles.end()));
        }

        deviceEnumeratedParticlesZipIter begin(){ return _begin;}
        deviceEnumeratedParticlesZipIter end() { return _end;}
    };

    typedef ParticleHostArray::iterator hostParticlesIter;
    typedef thrust::tuple<counter, hostParticlesIter> hostEnumeratedParticleIter;
    typedef thrust::zip_iterator<hostEnumeratedParticleIter> hostEnumeratedParticlesZipIter;

    class EnumeratedHostParticles
    {
        hostEnumeratedParticlesZipIter _begin;
        hostEnumeratedParticlesZipIter _end;
    public:
        EnumeratedHostParticles(ParticleHostArray& particles)
        {
            counter idx_begin(0);
            counter idx_end = idx_begin + particles.size();
            _begin = thrust::make_zip_iterator(thrust::make_tuple(idx_begin, particles.begin()));
            _end = thrust::make_zip_iterator(thrust::make_tuple(idx_end, particles.end()));
        }

        hostEnumeratedParticlesZipIter begin(){ return _begin;}
        hostEnumeratedParticlesZipIter end() { return _end;}
    };

    template <typename T>
    using enumeratedObject = thrust::tuple<uint,T>;

    template <typename T>
    class deviceArrayIter: public thrust::device_vector<T>::iterator{};

    template <typename T>
    class hostArrayIter: public thrust::host_vector<T>::iterator{};

    template <typename T>
    using enumeratedDeviceArrayIter = thrust::tuple<counter, deviceArrayIter<T>>;

    template <typename T>
    using enumeratedHostArrayIter = thrust::tuple<counter, hostArrayIter<T>>;

    template <typename T>
    using enumeratedDeviceArrayZipIter = thrust::zip_iterator<enumeratedDeviceArrayIter<T>>;

    template <typename T>
    using enumeratedHostArrayZipIter = thrust::zip_iterator<enumeratedHostArrayIter<T>>;

    template <typename T>
    class EnumeratedDeviceArray
    {
        enumeratedDeviceArrayZipIter<T> _begin;
        enumeratedDeviceArrayZipIter<T> _end;
    public:
        EnumeratedDeviceArray(thrust::device_vector<T>& arr)
        {
            counter idx_begin(0);
            counter idx_end = idx_begin + arr.size();
            _begin = thrust::make_zip_iterator(thrust::make_tuple(idx_begin, arr.begin()));
            _end = thrust::make_zip_iterator(thrust::make_tuple(idx_end, arr.end()));
        }

        enumeratedDeviceArrayZipIter<T> begin(){ return _begin;}
        enumeratedDeviceArrayZipIter<T> end() { return _end;}
    };
}

#endif