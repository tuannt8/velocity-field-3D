#ifndef VF_SEGMENT_MULTI_H
#define VF_SEGMENT_MULTI_H

#include "common/vf_common.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

namespace VF3D {
class Segment3DMulti{
    // The image
    int m_dim[3];
    std::vector<double> m_voxels; // row major

    std::vector<double> mean_intensity = {0, 0.36, 0.62, 0.73};
    double m_time_step = 5.;

    inline int index(int x, int y, int slice){
        return slice * m_dim[0] * m_dim[1] + x*m_dim[1] + y;
    }

    double interpolate_intensity(vec3 p){
        vec3i pti(std::floor(p[0]), std::floor(p[1]), std::floor(p[2]));
        vec3 relative_coord(p[0] - pti[0], p[1] - pti[1], p[2] - pti[2]);

        // TODO: write interpolate function
        double c00 = get_scalar(pti[0], pti[1], pti[2]) * (1 - relative_coord[0])
                    + get_scalar(pti[0] + 1, pti[1], pti[2]) * relative_coord[0];
        double c01 = get_scalar(pti[0], pti[1], pti[2]+1) * (1 - relative_coord[0])
                    + get_scalar(pti[0] + 1, pti[1], pti[2] + 1) * relative_coord[0];
        double c10 = get_scalar(pti[0], pti[1]+1, pti[2]) * (1 - relative_coord[0])
                    + get_scalar(pti[0] + 1, pti[1] + 1, pti[2]) * relative_coord[0];
        double c11 = get_scalar(pti[0], pti[1] + 1, pti[2]  +1) * (1 - relative_coord[0])
                    + get_scalar(pti[0] + 1, pti[1]+1, pti[2]+1) * relative_coord[0];


        double c0 = c00*(1-relative_coord[1]) + c10*relative_coord[1];
        double c1 = c01*(1-relative_coord[1]) + c11*relative_coord[1];

        double f =  c0*(1 - relative_coord[2]) + c1*relative_coord[2];
        return f;
    }

public:
    void init(std::string datapath){
        std::ifstream f(datapath, std::ios::in | std::ios::binary);

        // read dimension
        f.read((char*)&m_dim[0], sizeof(int));
        f.read((char*)&m_dim[1], sizeof(int));
        f.read((char*)&m_dim[2], sizeof(int));

        std::cout << "Loaded image size " << m_dim[0] << "x"
                  << m_dim[1] << "; "
                  << m_dim[2] << " slices" << std::endl;

        m_voxels.resize(m_dim[0]*m_dim[1]*m_dim[2]);
        std::vector<uint8_t> slice_data(m_dim[0]*m_dim[1]);
        double * ptr = m_voxels.data();
        for(int slice = 0; slice < m_dim[2]; slice++)
        {
            f.read((char*)slice_data.data(), slice_data.size());
            for(auto p : slice_data){
                *(ptr++) = (double)p / 255.;
            }
        }

        f.close();
    }

    vec3 get_displacement(vec3 pos, vec3 norm, int label_in, int label_out){
        assert(label_in >= 0 && label_in <4);
        assert(label_out >= 0 && label_out <4);

        double v = interpolate_intensity(pos);
        double c_in = mean_intensity[label_in];
        double c_out = mean_intensity[label_out];
        double f = (2*v - c_in - c_out) * (c_out - c_in);
        return -norm * f * m_time_step;
    }

    // Image info
    int no_slice(){return m_dim[2];}
    int im_width(){return m_dim[0];}
    int im_height(){return m_dim[1];}
    double get_scalar(int x, int y, int slice){
        // protect it
        if(x < 0 || x >= m_dim[0]
            || y < 0 || y >= m_dim[1]
            || slice < 0 || slice >= m_dim[2])
        {
            return mean_intensity[0];
        }

        return m_voxels[index(x, y, slice)];
    }

    // Row major
    const double * get_slice(int slice){
        return &m_voxels[slice * m_dim[0] * m_dim[1]];
    }
};
}
#endif // VF_SEGMENT_MULTI_H
