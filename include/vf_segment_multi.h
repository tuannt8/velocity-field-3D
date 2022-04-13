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
    double m_time_step = 0.5;
    double m_smooth_coef = 0.02;//relateable to intensity, avg = 0.02

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

    // Binary search for the exact location
    std::vector<vec3> get_displacement_triangle(vec3 pi0, vec3 pi1, vec3 pi2,
                                                vec3 norm, int label_in, int label_out, double max_distance)
    {
        std::vector<vec3> forces(3, vec3(0,0,0));

        vec3 p[3] = {pi0, pi1, pi2};
        // external force

        vec3 pos = (p[0] + p[1] + p[2])/3.0;

        // Binary search on normal direction
        vec3 p0 = pos - norm * max_distance;
        vec3 p1 = pos + norm * max_distance;

        for(int i = 0; i < 8; i++)
        {
            double f0 = get_priciple_external_force(p0, label_in, label_out);
            double f1 = get_priciple_external_force(p1, label_in, label_out);

            if(f0 > 0)// p0 is inside
            {
                if(f1 > 0){
                    // p1 is also inside
                    p0 = p1;
                    break;
                }
                else{
                    // p1 is outside
                    p1 = (p0 + p1)*0.5;
                }
            }
            else{ // p0 is outside
                p1 = p0; // move to p0
                break;
            }
        }

        // Final pos
        vec3 final_dis = (p0 + p1)*0.5 - pos;
        // If projection is small, switch to force model
        if(true
//                final_dis.norm() < max_distance/10
                ){
            return get_displacement_triangle(pi0, pi1, pi2, norm, label_in, label_out);
        }
        else
        {
            forces = {final_dis, final_dis, final_dis};
            return forces;
        }
    }

    std::vector<vec3> get_internal_force(vec3 pi0, vec3 pi1, vec3 pi2)
    {
        vec3 p[3] = {pi0, pi1, pi2};
        std::vector<vec3> forces(3, vec3(0,0,0));
        for(int i = 0; i < 3; i++){
            vec3 cp = p[i];
            vec3 others[2] = {p[(i+1)%3], p[(i+2)%3]};

            vec3 ab = others[1] - others[0];
            vec3 a = others[0];

            double t = ab.dot(a-cp) / std::pow(ab.norm(), 2);
            vec3 H = a + ab*t;
            vec3 HC = H - cp;

            vec3 f = -HC * ab.norm();

            forces[i] += f;
        }

        return forces;
    }

    std::vector<vec3> get_displacement_triangle(vec3 p0, vec3 p1, vec3 p2,
                                                vec3 norm, int label_in, int label_out)
    {
        vec3 p[3] = {p0, p1, p2};
        std::vector<vec3> forces(3, vec3(0,0,0));

        double area = std::pow( ((p[2] - p[0]).cross(p[1] - p[0])).norm(), 2 ) * 0.5;

        // external force
        std::vector<vec3> gauss_point = {vec3(0.3333, 0.33333, 0.33333)};
        std::vector<double> weight = {1.};
        for(size_t i = 0; i < gauss_point.size(); i++)
        {
            vec3 g = gauss_point[i];
            double w =  weight[i];
            vec3 pos = p[0]*g[0] + p[1]*g[1] + p[2]*g[2];
            vec3 f = get_displacement(pos, norm, label_in, label_out) * w * area;


            forces[0] += f * g[0];
            forces[1] += f * g[1];
            forces[2] += f * g[2];
        }


        // internal force
        for(int i = 0; i < 3; i++){
            vec3 cp = p[i];
            vec3 others[2] = {p[(i+1)%3], p[(i+2)%3]};

            vec3 ab = others[1] - others[0];
            vec3 a = others[0];

            double t = ab.dot(a-cp) / std::pow(ab.norm(), 2);
            vec3 H = a + ab*t;
            vec3 HC = H - cp;
            // HC.normalized();

            vec3 f = HC * ab.norm() * m_smooth_coef * m_time_step;

            forces[i] += f;
        }

        return forces;
    }

    double get_priciple_external_force(vec3 pos, int label_in, int label_out){
        assert(label_in >= 0 && label_in <4);
        assert(label_out >= 0 && label_out <4);

        double v = interpolate_intensity(pos);
        double c_in = mean_intensity[label_in];
        double c_out = mean_intensity[label_out];
        double f = -(2*v - c_in - c_out) * (c_out - c_in);

        return f;
    }

    vec3 get_external_force(vec3 pos, vec3 norm, int label_in, int label_out){
        assert(label_in >= 0 && label_in <4);
        assert(label_out >= 0 && label_out <4);

        double v = interpolate_intensity(pos);
        double c_in = mean_intensity[label_in];
        double c_out = mean_intensity[label_out];
        double f = (2*v - c_in - c_out) * (c_out - c_in);
        return -norm * f;
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
