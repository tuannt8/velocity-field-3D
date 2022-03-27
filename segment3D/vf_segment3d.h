#ifndef VF_SEGMENT3D_H
#define VF_SEGMENT3D_H

#include "common/vf_common.h"
#include <vector>
#include <string>
#include <fstream>

namespace VF3D {
class Segment3D{
    // The image
    int m_dim[3];
    std::vector<double> m_voxels; // row major

    double c_in = 0.83;
    double c_out = 0.39;
    double m_time_step = 2;//for final: 5 = max_dis = 1.5
    double m_smooth = 0.1;

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

    double get_curvature_tri(vec3 P, vec3 A, vec3 B)
    {
        // normal direction to l0
        // project l0 to {l1, l2}
        vec3 PA = A - P;
        vec3 PB = B - P;

        double cos_theta = PA.dot(PB) / (PA.norm() * PB.norm());

        double theta = acos(cos_theta);
        if(theta < 0){

        }

        return theta;
    }
    vec3 compute_curvature(vec3 pos, std::vector<std::vector<vec3>> tri_around){
        double c = 0;
        vec3 average_pos(0,0,0);
        for(auto poses : tri_around)
        {
            vec3 p1 = poses[1];
            vec3 p2 = poses[2];

            c += get_curvature_tri(pos, p1, p2);
            average_pos += (p1 + p2)*0.5;
        }
        average_pos /= tri_around.size();
        double curvature = std::max(0., 2*3.1415963 - c);

        vec3 curvature_vector = average_pos - pos;
        curvature_vector.normalize();

        return curvature_vector * curvature;
    }
public:
    Segment3D(){}
    ~Segment3D(){}

    void load(std::string path){
        std::ifstream f(path, std::ios::in | std::ios::binary);

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

    vec3 get_displacement(vec3 pos, vec3 norm,
                          std::vector<std::vector<vec3>> tris_around){
        // external force
        double v = interpolate_intensity(pos);
        double f = (2*v - c_in - c_out) * (c_out - c_in);
        vec3 f_ex = -norm * f;

        // internal force
        vec3 f_in = compute_curvature(pos, tris_around);

//        return f_in * m_smooth * m_time_step;

        return (f_ex + f_in * m_smooth)*m_time_step;
    }

    vec3 get_displacement(vec3 pos, vec3 norm){
        double v = interpolate_intensity(pos);
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
            return c_out;
        }

        return m_voxels[index(x, y, slice)];
    }

    // Row major
    const double * get_slice(int slice){
        return &m_voxels[slice * m_dim[0] * m_dim[1]];
    }
};
}



#endif // VF_SEGMENT3D_H
