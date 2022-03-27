#pragma once
#include <math.h>
#include <iostream>
#include "common/vf_common.h"

//
//  code adapted from El Topo's Enright driver:
//    https://github.com/tysonbrochu/eltopo/blob/master/talpa/drivers/enrightdriver.h
//

// 18.03.2022: Final release
// Nomore changes

namespace VF3D {
    class Enright{
        double m_time = 0;
        double m_dt = 0.01;

        void enright_velocity(double t, const vec3 pos, vec3 & out){
            double x = pos[0];
            double y = pos[1];
            double z = pos[2];

            out[0] = 2.0 * sin(    M_PI*x) * sin(    M_PI*x)
                    * sin(2.0*M_PI*y) * sin(2.0*M_PI*z);
            out[1] = -sin(2.0*M_PI*x) * sin(    M_PI*y)
                    * sin(    M_PI*y) * sin(2.0*M_PI*z);
            out[2] = -sin(2.0*M_PI*x) * sin(2.0*M_PI*y)
                    * sin(    M_PI*z) * sin(    M_PI*z);

            out *= sin(M_PI * t * 2 / 3);    // modulate with a period of 3
        }

    public:
        Enright(){

        };

        bool is_finish(){
            return m_time >= 3;
        }

        void next_step(){
            m_time += m_dt;
        }

        const vec3 get_displacement(const vec3 pos){
            vec3 v;
            enright_velocity(m_time, pos, v);
            vec3 k1 = v;
            enright_velocity(m_time + 0.5 * m_dt, pos + 0.5 * m_dt * k1, v);
            vec3 k2 = v;
            enright_velocity(m_time + 0.5 * m_dt, pos + 0.5 * m_dt * k2, v);
            vec3 k3 = v;
            enright_velocity(m_time + m_dt, pos + m_dt * k3, v);
            vec3 k4 = v;
            v = (1./6. * (k1 + k4) + 1./3. * (k2 + k3));

            return v*m_dt;
        }
    };
}


