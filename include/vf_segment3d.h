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
    double m_smooth = 0.02;

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


    std::vector<vec3> get_external_triangle_based(vec3 p0, vec3 p1, vec3 p2, vec3 norm)
    {
        static std::vector<vec3> gauss_point = {vec3(1./3., 1./3., 1./3.)};
        static std::vector<double> weight = {1.};

        std::vector<vec3> output(3, vec3(0,0,0));

        for(size_t i = 0; i < gauss_point.size(); i++)
        {
            vec3 g = gauss_point[i];
            double w = weight[i];

            vec3 pos = p0 * g[0] + p1 * g[1] + p2 * g[2];


            double v = interpolate_intensity(pos);
            double f = (2*v - c_in - c_out) * (c_out - c_in);
            vec3 f_v = -norm * (f * w);

            output[0] += f_v * g[0];
            output[1] += f_v * g[1];
            output[2] += f_v * g[2];
        }

        return output;
    }

    std::vector<vec3> get_internal_triangle_based(vec3 pi0, vec3 pi1, vec3 pi2)
    {
        vec3 p[3] = {pi0, pi1, pi2};
        std::vector<vec3> forces(3, vec3(0,0,0));
        for(int i = 0; i < 3; i++){
            vec3 cp = p[i];
            vec3 others[2] = {p[(i+1)%3], p[(i+2)%3]};

            vec3 ab = others[1] - others[0];
            vec3 a = others[0];

            double t = ab.dot(cp - a) / std::pow(ab.norm(), 2);
            vec3 H = a + ab*t;
            vec3 CH = H - cp;
            CH.normalize();

            vec3 f = CH * ab.norm();

            forces[i] = f * m_smooth;
        }

        return forces;
    }
public:
    Segment3D(){}
    ~Segment3D(){}

    bool is_finish(){
        static int iter = 0;
        return iter++ > 60;
    }

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

    std::vector<vec3> get_displacement_triangle_based(vec3 p0, vec3 p1, vec3 p2, vec3 norm)
    {
        auto external = get_external_triangle_based(p0, p1, p2, norm);
        auto internal = get_internal_triangle_based(p0, p1, p2);

        std::vector<vec3> out(3, vec3(0,0,0));

        for(int i = 0; i < 3; i++)
        {
            out[i] = (internal[i] + external[i])*m_time_step;
//            out[i] = ( internal[i])*m_time_step;
        }

        return out;
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
