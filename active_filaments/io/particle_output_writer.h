#ifndef __AH_SAVING__H__
#define __AH_SAVING__H__

#include <string>
#include <fstream>
#include <sstream>

#include "../particle.cuh"
#include "logging.h"

using namespace std;

namespace af 
{
    class ParticleOutputWriter
    {
    protected:
        uint print_count;
        ofstream* out;

        ParticleOutputWriter()
            : print_count(0), out(NULL){}

    public:

        ParticleOutputWriter(ofstream& stream)
            : print_count(0), out(&stream) {}

        virtual void update(const ParticleHostArray& particles) = 0;
    };
    

    class ParticleXYZWriter : public ParticleOutputWriter
    {
        const char types[10] = {'A','B','C','D','E','F','G','H','I','J'};

        const uint ppf;
        const uint prate;
        uint call_count;
        string file_name;
    public:
        ParticleXYZWriter(const string& filebase, uint print_rate = 1, uint prints_per_file = 0)
            : ppf(prints_per_file), prate(print_rate), 
            call_count(0)
        {
            file_name = filebase + ".xyz";
            out = new ofstream(file_name, ofstream::out);
            if (!out->is_open())
                throw exception("Could not open outputfile");
        }

        ~ParticleXYZWriter()
        {
            if (out->is_open()) out->close();
        }

        void update(const ParticleHostArray& particles)
        {
            // only print at the specified rate
            if (++call_count % prate != 1) return;

            // write output to a string stream
            stringstream xyz;
            xyz << particles.size() << endl;
            xyz << ++print_count << endl;
            for (auto p : particles)            
                xyz << "A " << p.r.x << " " << p.r.y << " " << p.r.z << endl;

            // write string stream
            (*out) << xyz.rdbuf();

            // TODO: Debug logger
            LOG(INFO) << "[ " << print_count << " ] Wrote particles to " << file_name;
        }
    };
}

#endif