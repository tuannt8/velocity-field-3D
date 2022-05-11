#ifndef VF_CURVATURE_H
#define VF_CURVATURE_H
#include "common/vf_common.h"

namespace VF3D{

enum speed_mode{
    speed_2x_edge_length = 0,
    speed_1x_edge_length = 1,
    speed_0_5x_edge_length = 2,
    speed_0_25x_edge_length = 3,
    speed_0_125x_edge_length = 4,
};

class Curvature_flow{
    double m_edge_length;
    double speed;
    double original_speed;
    int m_iter = 0;
    bool bDone = false;
public:
    void init(double avg_edge_legnth, speed_mode mode)
    {
        m_edge_length = avg_edge_legnth;
        std::map<speed_mode, double> speed_ratio = {
            {speed_2x_edge_length, 2},
            {speed_1x_edge_length, 1},
            {speed_0_5x_edge_length, 0.5},
            {speed_0_25x_edge_length, 0.25},
            {speed_0_125x_edge_length, 0.125}
        };
        speed =  speed_ratio[mode] * m_edge_length;
        original_speed = speed;
    }

    bool next(){
        m_iter++;

        if( m_iter * original_speed > 0.14)
        {
            speed = 0.14 - (m_iter - 1) * original_speed;
            std::cout << "===== Last speed " << speed << std::endl;

            if(speed < 0)
                return false;
        }

        return true;
    }


    vec3 get_new_position(vec3 pos, vec3 norm){

        vec3 newPos = pos + norm * speed;

        for(int i = 0; i < 3; i++)
        {
            newPos[i] = std::min(newPos[i], 1.);
            newPos[i] = std::max(newPos[i], 0.);
        }
        return newPos;
    }

};
}

#endif // VF_CURVATURE_H
