#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"
#include "json.hpp"
#include "spline.h"

#include "blackbox.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

Eigen::VectorXd optimize_JMT(double start_p, double start_v, double start_a, double end_p, double end_v, double end_a, double t) {
	
  Eigen::MatrixXd A = Eigen::MatrixXd(3,3);
  Eigen::VectorXd b = Eigen::VectorXd(3);
  Eigen::VectorXd x = Eigen::VectorXd(3);
  Eigen::VectorXd c = Eigen::VectorXd(6);

  double t2 = t * t;
  double t3 = t * t2;
  double t4 = t * t3;
  double t5 = t * t4;

    A << t3, t4, t5, 3*t2, 4*t3,  5*t4, 6*t , 12*t2, 20*t3;

    b << end_p - (start_p + start_v * t + 0.5 * start_a * t2),
         end_v - (start_v + start_a * t), end_a -  start_a;
		 
	Eigen::MatrixXd A_inverse = A.inverse();
	
	x = A_inverse * b;

    c << start_p, start_v, start_a * 0.5, x[0], x[1], x[2];
	
	std::cout << "coeffs:" << start_p << " " << start_v << " " << start_a * 0.5 << " " << x[0] << " " << x[1] << " " << x[2] << std::endl;
	std::cout << std::endl;
		
	return c;
}


double compute_JMT(Eigen::VectorXd c, double t) {

  double t0 = 1.0;
  double t2 = t * t;
  double t3 = t * t2;
  double t4 = t * t3;
  double t5 = t * t4;
  
  Eigen::VectorXd T = Eigen::VectorXd(6);
  T << t0, t, t2, t3, t4, t5;
  
  return T.transpose() * c;
}

double infinite_loop(double s, double limit) {
	
	if (s > limit / 2.0) {
		return std::fmod(s, limit);
		} else {return s;}
	
	}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

vector<double> d_lane(double x_value, double y_value, double yaw, vector<double> maps_x, vector<double> maps_y, vector<double> maps_s, double lane)
{
	vector<double> first_conversion = getFrenet(x_value, y_value, deg2rad(yaw), maps_x, maps_y);
	double s_value = first_conversion[0];
	double d_value = first_conversion[1] + (2.0 + lane * 4.0);
	vector<double> coords = getXY(s_value, d_value, maps_s, maps_x, maps_y);
	//std::cout << "INPUT: " << x_value << ", " << y_value << ", " << yaw << " @ " << lane << std::endl;
	//std::cout << "OUTPUT: " << s_value << ", " << d_value <<  std::endl;
	return {coords[0]-x_value, coords[1]-y_value};
	}
	
	

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
  
  // Integrating waypoints data with continous information at start and end 
  // This will help rendering a perfect loop
	for (int i = (map_waypoints_s.size()-1) ; i < (map_waypoints_s.size()-50); i--) { 
	  float s_2 = map_waypoints_s[i] - max_s;
	  double x_2 = map_waypoints_x[i];
	  double y_2 = map_waypoints_y[i];
	  float d_x_2 = map_waypoints_dx[i];
	  float d_y_2 = map_waypoints_dy[i];
	  map_waypoints_s.insert(map_waypoints_s.begin(), s_2);
	  map_waypoints_x.insert(map_waypoints_x.begin(), x_2);
	  map_waypoints_y.insert(map_waypoints_y.begin(), y_2);
	  map_waypoints_dx.insert(map_waypoints_dx.begin(), d_x_2);
	  map_waypoints_dy.insert(map_waypoints_dy.begin(), d_y_2);
	}

	for (int i = 0; i < 50; i++) { 
	  float s_2 = map_waypoints_s[i] + max_s;
	  double x_2 = map_waypoints_x[i];
	  double y_2 = map_waypoints_y[i];
	  float d_x_2 = map_waypoints_dx[i];
	  float d_y_2 = map_waypoints_dy[i];
	  map_waypoints_s.push_back(s_2);
	  map_waypoints_x.push_back(x_2);
	  map_waypoints_y.push_back(y_2);
	  map_waypoints_dx.push_back(d_x_2);
	  map_waypoints_dy.push_back(d_y_2);
	}
  
  
  // Blackbox is a helper class to help us record and store what is useful 
  Blackbox trip;

  // General parameters for the simulation
  trip.speed_limit = 45;
  trip.preferred_lane = 1.0;
  trip.running_lane = 1.0;
  trip.target_lane = 1.0;
  trip.target_speed = trip.speed_limit;
  trip.last_s = -1;
  trip.last_d = -1;
  trip.cold_start = 1;
  
  // Modelling the path of the circuit using spline functions
  trip.spline_x.set_points(map_waypoints_s, map_waypoints_x);      // fit a spline for x coordinates based on s
  trip.spline_y.set_points(map_waypoints_s, map_waypoints_y);      // fit a spline for y coordinates based on s
  trip.spline_dx.set_points(map_waypoints_s, map_waypoints_dx);    // fit a spline for dx coordinates based on s
  trip.spline_dy.set_points(map_waypoints_s, map_waypoints_dy);    // fit a spline for dy coordinates based on s
  
  // Information on the waypoints
  trip.max_s = max_s;
  
  // There are 181 waypoints
  trip.lenght_waypoints = map_waypoints_x.size();
  
  h.onMessage([&trip,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
		
		std::mt19937 rng;
		rng.seed(std::random_device()());
		std::uniform_int_distribution<std::mt19937::result_type> switch_lane(0,9); // distribution in range [0, 9]
		std::uniform_int_distribution<std::mt19937::result_type> which_lane(0,2); // distribution in range [0, 2]
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
			
			// Cold start
			if (trip.last_s == -1) {trip.last_s = car_s;}
			if (trip.last_d == -1) {trip.last_d = car_d;}

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
			
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
			
			// Fast evaluation of the sensor data 
			trip.lane_front = {999, 999, 999};
			trip.lane_rear = {999, 999, 999}; 
			trip.crowd = {0, 0, 0};
			trip.front_speed = {trip.speed_limit, trip.speed_limit, trip.speed_limit};
			trip.back_speed = {0, 0, 0};
	
			for (int i = 0; i < sensor_fusion.size(); i++) {

				  int id = sensor_fusion[i][0];
				  double s = sensor_fusion[i][5];
				  double d = sensor_fusion[i][6];
				  double vx = sensor_fusion[i][3];
				  double vy = sensor_fusion[i][4];
				  			    
				  double speed = sqrt(vx * vx + vy * vy) * 2.0;
				  
				  double distance_from_us = abs(car_s - s);
				  double front =  ((s-car_s) > 0) - ((s-car_s) < 0);
				  double lane = -1.0;
				  
				  for (int j = 0; j < 3; j++) {
					  
					if ((d >= 0.0 + 4.0 * j) & (d < 4.0 + 4.0 * j)) {
					  lane = j;
					  if (front > 0) {
						  if (distance_from_us < trip.lane_front[j]) {
							  trip.lane_front[j] = distance_from_us;
							  trip.front_speed[j] = speed;
							  }
						  trip.crowd[j] += 1;
					  } else {
						  if (distance_from_us < trip.lane_rear[j]) {
							  trip.lane_rear[j] = distance_from_us;
							  trip.back_speed[j] = speed;
						  }
						  }
					  }
					  
				  }
				  
				  // std:cout << "car: " << id << " s " << s << " d " << d << " lane " << lane << " speed " << speed << " distance " << distance_from_us << std::endl;
				  }
				
				std::cout << "L > " << trip.lane_front[0] << " | " << trip.lane_rear[0] << " cars: " << trip.crowd[0] << " cruise speed: " << trip.front_speed[0] << std::endl;
				std::cout << "M > " << trip.lane_front[1] << " | " << trip.lane_rear[1] << " cars: " << trip.crowd[1] << " cruise speed: " << trip.front_speed[1] << std::endl;
				std::cout << "R > " << trip.lane_front[2] << " | " << trip.lane_rear[2] << " cars: " << trip.crowd[2] << " cruise speed: " << trip.front_speed[2] << std::endl;
			
			// Behaviour planner
			double opportunity_cost;
			double feasibility_cost;
			std::vector<double> total_cost = {0.0, 0.0, 0.0};
			
			double manoeuvre_space = 20;
			double security_distance = 30;
			
			// Updating running lane record based on last target reached
			trip.running_lane = trip.target_lane;
			double next_action = trip.running_lane;
			double best_action = 999.0;
			 
			for (int j = 0; j < 3; j++) {
				if ((trip.lane_rear[j] > manoeuvre_space * 0.75) & (trip.lane_front[j] > manoeuvre_space)) {
					feasibility_cost = 0.0;
					} else {feasibility_cost = 9.0;}
				
				opportunity_cost = trip.crowd[j] * 0.2 + 30.0 / trip.lane_front[j];
				
				total_cost[j] = opportunity_cost + feasibility_cost;
				
				//std::cout << "Feasibility " << feasibility_cost << std::endl;
				//std::cout << "Opportunity " << opportunity_cost << std::endl;
				
				std::cout << "Cost for running in lane " << j << " : " << total_cost[j] << std::endl;
			}
			
			std::cout << "Present speed: " << car_speed << std::endl;
			std::cout << "Running in lane " << trip.running_lane << std::endl;
			
			// Correction for double lane change
			if (trip.running_lane == 0) {total_cost[2] += total_cost[1];}
			if (trip.running_lane == 2) {total_cost[0] += total_cost[1];}
			
			// Correction for lane permanence
			total_cost[trip.running_lane] = total_cost[trip.running_lane] * 0.8;
			
			// Finding the best thing to do
			for (int j = 0; j < 3; j++) {
				if ((total_cost[j] < best_action) & (total_cost[j] < 1.5)) {
					next_action = j;
					best_action = total_cost[j];					
					}
			}
			
			if (next_action == trip.running_lane) {
				std::cout << "We keep the lane " << std::endl;
				if (trip.lane_front[next_action] < security_distance) {
					//trip.target_speed = max(0.0, car_speed - slow_down);
					std::cout << "SLOWER CAR IN RANGE " << std::endl;
					trip.target_speed = trip.front_speed[trip.running_lane];
					std::cout << "But we have to slow down at speed : " << trip.target_speed << " mph" << std::endl;
					} else {
						trip.target_speed = trip.speed_limit;
						std::cout << "At our cruise speed : " << trip.target_speed << " mph" << std::endl;
						}
				} else {
					std::cout << "We move to lane " << next_action << std::endl;
					trip.target_lane = next_action;
					}
			
			// Fixing the next target waypoint
			int waypoint = NextWaypoint(car_x, car_y, car_yaw, map_waypoints_x, map_waypoints_y);
			
			// check not to get struck exactly on a waypoint
			if (((car_x == map_waypoints_x[waypoint]) & (car_y == map_waypoints_y[waypoint])) | (waypoint==trip.lenght_waypoints)) {				
				if (waypoint==trip.lenght_waypoints) {
					waypoint = 1;
				} else {
					waypoint += 1;
					}
			}
						
			// Reporting th situation with the car
			//std::cout << "NEXT WAYPOINT no " << waypoint << " x: " << map_waypoints_x[waypoint] << " y: " << map_waypoints_y[waypoint] <<  std::endl;
			
			std::cout << "OVERVIEW" << std::endl;
			std::cout << "car x " << car_x << " car y " << car_y << std::endl;
			std::cout << "car s " << car_s << " car d " << car_d << std::endl;
			std::cout << "car yaw " << car_yaw << std::endl; 
			std::cout << "car actual speed " << car_speed << " car target speed " << trip.target_speed << std::endl;
			std::cout << "car actual lane "  << trip.running_lane << " car target lane " << trip.target_lane << std::endl;
 
          	json msgJson;

          	vector<double> next_x_vals, next_y_vals;
			vector<double> rough_s_vals, rough_x_vals, rough_y_vals;
			vector<double> smooth_s_vals, smooth_x_vals, smooth_y_vals;
			
			// Re-using the previous trajectory
			if (previous_path_x.size() > 0) {
				
				for(int i = 0; i < previous_path_x.size(); i++)
					{
						  next_x_vals.push_back(previous_path_x[i]);
						  next_y_vals.push_back(previous_path_y[i]);
					}
			}
			
          	// defining a path made up of (x,y) points that the car will visit sequentially every .02 seconds
			
			int limit_to_reuse = 15;
			if (previous_path_x.size() > limit_to_reuse) {
				
				std::cout << "Reusing previous trajectory..." << std::endl;
				
				} else {
					
					std::cout << "residuals: " << previous_path_x.size() << std::endl;
					
					// Computing a new trajectory
					// Configuring up JMT
					
					double interval = 0.02;
					double time_involved;
					double target_s;
					double target_d;
					double target_speed_ms;
					double number_points;
					double car_speed_ms;
					double average_speed_ms;
					
					if (trip.cold_start==1) {
						
						std::cout << "COLD START" << std::endl;
						
						target_d = car_d;
						time_involved = 2.25;
						target_speed_ms = trip.target_speed * 1000 / 60 / 60 * 1.60934 * 0.5;
						target_s = car_s + (target_speed_ms * time_involved) * 0.5;
						number_points = time_involved / interval;
						
						trip.cold_start = 0;
						
					} else {
					
						time_involved = 2.25;
						
						if (trip.target_lane == 0) {
								target_d = 2.25;
								} else {
									if (trip.target_lane==1){
										target_d = 6;
										} else {
											target_d = 9.75;
											} 
								}
						
						target_speed_ms = trip.target_speed * 1000 / 60 / 60 * 1.60934;
						car_speed_ms = car_speed * 1000 / 60 / 60 * 1.60934;
						average_speed_ms = (car_speed + trip.target_speed) * 0.5 * 1000 / 60 / 60 * 1.60934;
						target_s = trip.last_s + average_speed_ms * time_involved;
						
						number_points = time_involved / interval;
						
					}
					
					std::cout << "start s " << trip.last_s << " end s " << target_s << std::endl;
					std::cout << "start d " << trip.last_d << " end d " << target_d << std::endl;
					std::cout << "actual speed " << car_speed_ms << " target speed " << target_speed_ms << std::endl;
					
					Eigen::VectorXd JMT_s = optimize_JMT(trip.last_s, car_speed_ms, 0.0, target_s, target_speed_ms, 0.0, time_involved);
					Eigen::VectorXd JMT_d = optimize_JMT(trip.last_d, 0.0, 0.0, target_d, 0.0, 0.0, time_involved);
					
					double s1 = compute_JMT(JMT_s, 0.0);
					double s2 = compute_JMT(JMT_s, time_involved);
					double d1 = compute_JMT(JMT_d, 0.0);
					double d2 = compute_JMT(JMT_d, time_involved);
					
					std::cout << "s FROM " << s1 << " TO " << s2 << std::endl;
					std::cout << "d FROM " << d1 << " TO " << d2 << std::endl;
					
					// Converting JMT fitted polynomial into path of global coordinates
					std::cout << " " << std::endl;
					
					for (double t = 1; t < number_points; t++) {       // Please note t = 1 because of latency
							
							double new_s = infinite_loop(compute_JMT(JMT_s, t * interval), trip.max_s);
							trip.last_s = new_s;
							double new_d = compute_JMT(JMT_d, t * interval);
							trip.last_d = new_d;
							
							double x_edge = trip.spline_x(new_s);
							double y_edge = trip.spline_y(new_s);
							double dx = trip.spline_dx(new_s);
							double dy = trip.spline_dy(new_s);
							
							double new_x = x_edge + dx * new_d;
							double new_y = y_edge + dy * new_d;
														
							next_x_vals.push_back(new_x);
							next_y_vals.push_back(new_y);
							
							//std::cout << t * interval << " s " << new_s << " d " << new_d << " x " << new_x << " y " << new_y << std::endl;
						}
				}
				
			std::cout << " " << std::endl;
			
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}