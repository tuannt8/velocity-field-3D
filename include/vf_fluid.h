#ifndef FLUID_VTK_H
#define FLUID_VTK_H

#include <string>
#include "common/vf_common.h"
#include <fstream>

#include "common/kdtree/KDTreeVectorOfVectorsAdaptor.h"
#include "common/Eigen/Eigen"

// Load from VTK file
namespace VF3D {
typedef KDTreeVectorOfVectorsAdaptor< vec3array, double, 3,
nanoflann::metric_L2_Simple>  kdtree_flann;

class Anisotropic{
    // influence radius
    float m_h = 0.025;//0.03;
    float m_r = 0.05;//0.05;
    float sigma = 1;
    std::vector<vec3> m_particlePos;

    // It is final now
    double kn = 0.5;//20.;//0.5; // The smaller, the smaller blob
    double kr = 6.;//6.; // the larger, the more elipsoid
    double ks = 4000.;// 4000.; // the smaller ks, the tighter? 4000 seem right
    // Optimal 6000
    // ks / kr should be
    // 0.001 * ks / kr -> kn
    size_t no_nei_thres = 15; // No need for anisotropic

    std::vector<mat3> m_G; // anisotropic matrix
    std::vector<double> m_detG;
    std::vector<bool> m_bKernelBuilt;

    std::shared_ptr<kdtree_flann> m_kdtree;
    void build_kd_tree(){
        m_kdtree = std::shared_ptr<kdtree_flann>(
                    new kdtree_flann(3,
                                     m_particlePos,
                                     10)); // leaf size
        m_kdtree->index->buildIndex();
    }

    std::vector<std::pair<size_t, double>> neighbor_search(vec3 pos,
                                                           double distance){
        std::vector<std::pair<size_t, double>> resultSet;
        double distance2 = distance * distance;
        m_kdtree->index->radiusSearch(pos.data(),
                                      distance2,
                                      resultSet,
                                      nanoflann::SearchParams());

        return resultSet;
    }

    void laplace_smooth(){
        double lamda = 0.9;// 0.7; // closer to 1 -> more smooth
        std::vector<vec3> smoothed_particles = m_particlePos;

        for (size_t i = 0; i < m_particlePos.size(); i++)
        {
            auto curPos = m_particlePos[i];

            // Search for neighbor particles
            auto neighbors = neighbor_search(curPos, m_r);

            if (neighbors.size() > 0)
            {
                vec3 new_pos(0, 0, 0);
                double sum_omega = 0;
                for (auto nidx : neighbors)
                {
                    auto & pj = m_particlePos[nidx.first];
                    double omega = weight_func(curPos, pj, m_r);
                    sum_omega += omega;
                    new_pos += pj*omega;
                }
                new_pos /= sum_omega;
                new_pos = new_pos*lamda + curPos*(1-lamda);

                smoothed_particles[i] = new_pos;
            }
        }
        m_particlePos = smoothed_particles;
    }

    double weight_func(vec3 p1, vec3 p2, double radius){
        double r = (p1 - p2).norm() / radius;
        if (r > 1)
        {
            return 0;
        }

        return 1 - r*r*r;
    }

    const mat3 &get_transform_mat(int idx){
        if(!m_bKernelBuilt[idx])
        {
            compute_tranformation_mat_for_particle(idx);
            m_bKernelBuilt[idx] = true;
        }

        return m_G[idx];
    }

    void compute_tranformation_mat_for_particle(int idx)
    {
        // 1 - 6 - 6000 - 10
        vec3 pos = m_particlePos[idx];
        auto neighbor = neighbor_search(pos, m_r);

        // connected only, i.e. r < m_h (0.0195)
        std::vector<std::pair<size_t, double>> real_neighbor;
        for(auto n : neighbor)
        {
            if(m_labels[n.first] == m_labels[idx])
                real_neighbor.push_back(n);
        }
        neighbor = real_neighbor;

        // Weighted mean position
        vec3 x_w(0., 0., 0.);
        double sum_omega = 0.;
        for(auto ni : neighbor)
        {
            vec3 nPos = m_particlePos[ni.first];
            double omega = weight_func(pos, nPos, m_r);
            x_w += nPos * omega;
            sum_omega += omega;
        }
        x_w /= sum_omega;

        // Orientation matrix
        mat3 C;
        C.setZero();
        for(auto ni : neighbor)
        {
            auto nPos = m_particlePos[ni.first];
            double omega = weight_func(x_w, nPos, m_r);
            vec3 r_w = nPos - x_w;

            mat3 oo = r_w * r_w.transpose();
            C = C + oo * omega;
        }
        C /= sum_omega;

        // SVD solving
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(C, Eigen::ComputeFullU);
        mat3 Q = svd.matrixU();
        vec3 L = svd.singularValues();

        mat3 Sigma = mat3::Zero();

        for (int d = 0; d < 3; d++)
        {
            if(neighbor.size() > no_nei_thres)
            {
                Sigma(d,d) = std::max(L[d], L[0] / kr) * ks;
            }
            else
            {
                Sigma(d,d) = 1.0 * kn;
            }
        }

        mat3 Sigma_inv = mat3::Zero();
        Sigma_inv(0, 0) = 1.0 / Sigma(0, 0);
        Sigma_inv(1, 1) = 1.0 / Sigma(1, 1);
        Sigma_inv(2, 2) = 1.0 / Sigma(2, 2);

        m_G[idx] = (1./(m_h)) *  Q * Sigma_inv * Q.transpose();
        m_detG[idx] = m_G[idx].determinant();
    }
public:
    std::vector<uint8_t> m_labels;
    void build_labeling(){

    }
public:
    double scalar_value(vec3 pos, int idx){
        auto G = get_transform_mat(idx);
        vec3 aniso_radius = G * (pos - m_particlePos[idx]);
        double rah = aniso_radius.norm();

        if(rah > 1)
            return 0;
        else
            return sigma
                * 3.3e-6
                * m_detG[idx]
                * (1 - rah*rah*rah);
    }
    double scalar_value(vec3 pos){
        auto neighbor = neighbor_search(pos, m_r);

        double out = 0;
        for(auto ni : neighbor)
        {
            float sv = scalar_value(pos, ni.first);
            out += sv;
        }

        return out;
    }

public:
    Anisotropic(){

    };
    ~Anisotropic(){}

public: // Debugging
    void build( std::vector<vec3> particlePos, std::vector<uint8_t> labels)
    {
        m_particlePos = particlePos;
        m_labels = labels;
        build_kd_tree();

        // Allocate
        m_G.resize(m_particlePos.size());
        m_detG.resize(m_particlePos.size());

        // Reset kernel
        m_bKernelBuilt = std::vector<bool>(m_particlePos.size(), false);
    }
};

class fluid
{
    std::shared_ptr<std::ifstream>  m_file_stream;
    int num_particle;
    double spacing_delta;
    int m_speed = 1;


    double m_r = 0.039; // 1.3 * 0.0195

    std::shared_ptr<kdtree_flann> m_kdtree;
    std::vector<vec3> m_curParticlesPos;
    std::vector<vec3> m_nextParticlesPos;


    double weight_func(vec3 pos1, vec3 pos2, double radius){
        double r = (pos1 - pos2).norm() / radius;
        if (r > 1){
            return 0;
        }

        return 1 - r*r*r;
    }

    std::vector<std::pair<size_t, double>> find_neighbor(vec3 pos){
        std::vector<std::pair<size_t, double>> resultSet;
        m_kdtree->index->radiusSearch(pos.data(),
                                      m_r * m_r,
                                      resultSet,
                                      nanoflann::SearchParams());
        return resultSet;
    }


    void clamp_pos(vec3 &p)
    {
        for(int i = 0; i < 3; i++)
        {
            p[i] = std::max(0., p[i]);
            p[i] = std::min(1., p[i]);
        }
    }

    // Anisotropic kernel
    Anisotropic m_kernel;
    vec3 get_projection_displacement(vec3 pos, vec3 norm){
        double alpha = 0.001;

        double phi = m_kernel.scalar_value(pos);

        if(phi < 1e-8){ // outside
            return -norm * 0.015 * 0.1;
        }
        else{
            return norm * phi * alpha;
        }
    }

public:
    fluid(){};
    ~fluid(){}

    // Maximum displacement between two iterations is 0.008
    // By loading multiple iteration at one, the maximum displacement
    //      is multiplied and increases
    void set_speed(int speed){
        m_speed = speed;
        std::cout <<"Fluid speed " << speed << std::endl;
    }

    vec3 get_displacement(vec3 pos, vec3 norm){
        auto neighbors = find_neighbor(pos);

        if(neighbors.size() == 0)
        {
            return -norm * m_r / 10;
        }

        vec3 displace(0,0,0);
        double sum_omega = 0;
        for(auto ni : neighbors)
        {
            auto npos = m_curParticlesPos[ni.first];
            double omega = weight_func(pos, npos, m_r);

            sum_omega += omega;
            displace += (m_nextParticlesPos[ni.first]
                    - m_curParticlesPos[ni.first])*omega;
        }

        assert(sum_omega > 1e-10);
        displace /= sum_omega;

        vec3 newPos = pos + displace;

        // Projection
        vec3 project = get_projection_displacement(newPos, norm);
        newPos = newPos + project;

        clamp_pos(newPos);

        displace = newPos - pos;

        return displace;
    }



    void init(std::string path){
        m_file_stream = std::shared_ptr<std::ifstream>(
                    new std::ifstream(path, std::ios::in | std::ios::binary));
        if(!m_file_stream->is_open()){
            std::cout << "Cannot open data file: " << path;
            exit(1);
        }

       m_file_stream->read((char*)&num_particle, sizeof(int));
       m_file_stream->read((char*)&spacing_delta, sizeof(double));

       std::cout << num_particle << " * " << spacing_delta << std::endl << std::flush;

       // check header num_particle and spacing delta
       if(num_particle != 79036){
           std::cout << "Wrong file data for fluid 3D";
           exit(1);
       }


    }

    bool next_iteration(){
        if(m_nextParticlesPos.size() == 0)
            m_curParticlesPos = load_next_particle();
        else
            m_curParticlesPos = m_nextParticlesPos;

        m_nextParticlesPos = load_next_particle(m_speed);

        // Build annisotropic kernel
        m_kernel.build(m_nextParticlesPos, m_labels);

        // Build kd-tree
        m_kdtree = std::shared_ptr<kdtree_flann>(
                    new kdtree_flann(3,
                                     m_curParticlesPos,
                                     10)); // leaf size
        m_kdtree->index->buildIndex();

        return m_nextParticlesPos.size() > 0;
    }

    float get_cur_time(){return m_cur_time;}
    size_t get_fluid_sim_step(){return m_cur_time/0.0005;}

    // Offset in data file
    void offset(int num_iter){
        size_t buffer_size = sizeof(float) // timestep
                + num_particle * 3 * sizeof(unsigned short) // position
                + num_particle; // buffer connected component
        m_file_stream->seekg(buffer_size*(num_iter), std::ios_base::cur);
    }
private:
    float m_cur_time;
    std::vector<uint8_t> m_labels;
    std::vector<vec3> load_next_particle(int speed = 1){
        // Offset data
        offset(speed - 1);

        std::vector<short> pos_i(num_particle * 3);
        m_labels.resize(num_particle);

        if(!m_file_stream->read((char*)&m_cur_time, sizeof(float)) ||
           !m_file_stream->read((char*)pos_i.data(),pos_i.size() * sizeof(unsigned short)) ||
           !m_file_stream->read((char*)m_labels.data(), num_particle))
            return std::vector<VF3D::vec3>();

        vec3array pos(num_particle);
        for(int i = 0; i < num_particle; i++){
            pos[i] = vec3((unsigned short)pos_i[i*3],
                    (unsigned short)pos_i[i*3 + 1],
                    (unsigned short)pos_i[i*3 + 2]) * spacing_delta;
        }

        return pos;
    }
};
}
#endif // FLUID_VTK_H
