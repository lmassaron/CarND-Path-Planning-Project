#ifndef BLACKBOX_H
#define BLACKBOX_H
#include <vector>
#include "spline.h"

class Blackbox {
public:
  
  /*
  * Variables
  */
  
  double speed_limit;
  double target_speed;
  double preferred_lane;
  double running_lane;
  double target_lane;
  double last_s;
  double last_d;
  double max_s;
  int lenght_waypoints;
  int cold_start;
  std::vector<double> lane_front;
  std::vector<double> lane_rear; 
  std::vector<double> crowd;
  std::vector<double> front_speed;
  std::vector<double> back_speed;
  tk::spline spline_x, spline_y, spline_dx, spline_dy;
  
  /*
  * Initialize The class.
  */
  void Init();
  
};

#endif /* BLACKBOX_H */