#ifndef FLUID_VTK_H
#define FLUID_VTK_H

#include <string>
#include "common/vf_common.h"
#include "vf_anisotropic.h"
#include <fstream>

// Load from VTK file
namespace VF3D {
class fluid
{
    std::ifstream * m_file_stream = nullptr;
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
    ~fluid(){release();}

    void set_speed(int speed){m_speed = speed;}


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


    // New data type
    void init(std::string path){
        m_file_stream = new std::ifstream(path, std::ios::in | std::ios::binary);
        if(!m_file_stream->is_open()){
            throw "Cannot open data file";
        }

       m_file_stream->read((char*)&num_particle, sizeof(int));
       m_file_stream->read((char*)&spacing_delta, sizeof(double));

       std::cout << num_particle << " * " << spacing_delta << std::endl << std::flush;
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
    void release(){
        if(m_file_stream){
            m_file_stream->close();
            delete m_file_stream;
            m_file_stream = nullptr;
        }
    }

    void load_anisotropic(int iter){
        size_t buffer_size = sizeof(float) // timestep
                + num_particle * 3 * sizeof(unsigned short)
                + num_particle; // buffer
        size_t begin_size = sizeof(int) + sizeof(double);

        m_file_stream->seekg(buffer_size * (iter) + begin_size);

        auto pos = load_next_particle();
        std::cout << "Load vel at iter " << iter << ", " << m_cur_time << std::endl;

        return;

        m_kernel.build(pos, m_labels);

//        // Log it
//        {
//            std::stringstream name;
//            name << "LOG/particle_" << iter << "_"
//                 << (int)(get_fluid_sim_step()) << ".obj";

//            std::cout << "Log particles " << name.str() << std::endl;

//            std::ofstream f(name.str());
//            for(auto p : pos)
//                f << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
//            f.close();
//        }
    }

//    std::vector<short> pos_i;


    float get_cur_time(){return m_cur_time;}
    size_t get_fluid_sim_step(){return m_cur_time/0.0005;}

private:
    float m_cur_time;
    std::vector<uint8_t> m_labels;
    std::vector<vec3> load_next_particle(int speed = 1){
        size_t buffer_size = sizeof(float) // timestep
                + num_particle * 3 * sizeof(unsigned short) // position
                + num_particle; // buffer connected component
        m_file_stream->seekg(buffer_size*(speed-1), std::ios_base::cur);

        std::vector<short> pos_i(num_particle * 3);
        m_labels.resize(num_particle);

        if(!m_file_stream->read((char*)&m_cur_time, sizeof(float)) ||
           !m_file_stream->read((char*)pos_i.data(),pos_i.size() * sizeof(unsigned short)) ||
           !m_file_stream->read((char*)m_labels.data(), num_particle))
            return std::vector<VF3D::vec3>();

//        std::cout << m_cur_time << std::endl << std::flush;

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
