/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#define PI 3.14159265

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  if(is_initialized)
  {
   return ; 
  }

  // Setting the number of particles - Tuning for speed vs accuracy
  num_particles = 500;

  default_random_engine gen; // Random number gen

  // This Line creates a Normal (Gaussian) Distribition for x,y,theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for(int i=0; i <num_particles;i++)
  {
   Particle particle;
   particle.id = i;
   particle.x = dist_x(gen);
   particle.y = dist_y(gen);
   particle.theta = dist_theta(gen);
   particle.weight = 1;

   particles.push_back(particle);
  }

  is_initialized = true;// Initialization done

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  default_random_engine gen;
  
  for(int i=0;i<num_particles;i++)
  {
     
   // CTRV Motion model update equations for position x, position y and theta
   if(fabs(yaw_rate)<0.001f)
   {
    particles[i].x += velocity*delta_t*cos(particles[i].theta);
    particles[i].y += velocity*delta_t*sin(particles[i].theta);
   }
   else
   {
    particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) -sin(particles[i].theta));
    particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta) -cos(particles[i].theta + yaw_rate*delta_t));
   }
   particles[i].theta +=  yaw_rate*delta_t;

   // Adding some Gaussian noise
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  particles[i].x += dist_x(gen);
  particles[i].y += dist_y(gen);
  particles[i].theta += dist_theta(gen);
  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

 unsigned int num_of_obs = observations.size();
 unsigned int num_of_predictions = predicted.size();
 
 for(int i=0;i<num_of_obs;i++)
 {
  // For each obserbation we need to find the least distance to landmark and associate obs with Landmark id
  double min_dist = numeric_limits<double>::max();
  
  // Map Id init
  int mapId = -1;// Basically map id not found in this case

  for(int j=0;j<num_of_predictions;j++)
  {
   double dist_x = observations[i].x - predicted[j].x;
   double dist_y = observations[i].y - predicted[j].y;
   
   double distance = dist_x*dist_x + dist_y*dist_y;
   
   if(distance < min_dist)
   {
    min_dist = distance;
    mapId = predicted[j].id;
   }
  }

  // Only way to return this to the other function will be to update the id in the class
  observations[i].id = mapId; 
 }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


 double Lndmrk_range = std_landmark[0];
 double Lndmrk_phi   = std_landmark[1];
 
 for(int i=0;i<num_particles;i++)
 {
  
  double x = particles[i].x;
  double y = particles[i].y;
  double theta = particles[i].theta;

  //Find landmarks in particle range
  vector<LandmarkObs> LandmrksInRange;
  for(int j=0;j<map_landmarks.landmark_list.size();j++)
  {
   double LndmrkX = map_landmarks.landmark_list[j].x_f;
   double LndmrkY = map_landmarks.landmark_list[j].y_f;
   int id = map_landmarks.landmark_list[j].id_i;
   
   double del_x = x - LndmrkX;
   double del_y = y - LndmrkY;
 
   if(del_x*del_x + del_y*del_y <= sensor_range*sensor_range)
   {
    LandmrksInRange.push_back(LandmarkObs{id,LndmrkX,LndmrkY});
   }

  }

  // Transform observation coordinates to Map coordinates
  vector<LandmarkObs> ObsMap;
  for(int k=0; k <observations.size(); k++)
  {
   double x_conv = cos(theta)*observations[k].x - sin(theta)*observations[k].y + x;
   double y_conv = sin(theta)*observations[k].x + cos(theta)*observations[k].y + y;
   ObsMap.push_back(LandmarkObs{observations[k].id,x_conv,y_conv});
  }

  // Data association to Landmark
  dataAssociation(LandmrksInRange,ObsMap);
  
  // Weights - Reset and caluclate new weights
  particles[i].weight = 1.0f;

  //Calculate New Weights
  for(int m=0;m < ObsMap.size();m++)
  {
    double ObsX = ObsMap[m].x;
    double ObsY = ObsMap[m].y;
    int landmarkid = ObsMap[m].id;

    double assoc_x,assoc_y;
    bool found =false;
    int k=0;
    while (!found && k <LandmrksInRange.size())
    {
     if(LandmrksInRange[k].id == landmarkid)
     {
      found = true;
      assoc_x = LandmrksInRange[k].x;
      assoc_y = LandmrksInRange[k].y;
     }
     k++;
    }

   // Apply the gaussian distribution
   double x_diff = ObsX - assoc_x;
   double y_diff = ObsY - assoc_y;
   double weight = (1/(2*PI*Lndmrk_range*Lndmrk_phi))*exp(-(x_diff*x_diff/(2*Lndmrk_range*Lndmrk_range)) - (y_diff*y_diff/(2*Lndmrk_phi*Lndmrk_phi)));
   
   particles[i].weight *= weight;

  }
 }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

default_random_engine gen;

// Get Max weight
vector<double> weights;
double max_weight = numeric_limits<double>::min();

for(int i=0;i<num_particles;i++)
{
 weights.push_back(particles[i].weight);
 if(particles[i].weight>max_weight)
 {
 max_weight = particles[i].weight; 
 }
}

// Creating the distribution
uniform_real_distribution<double> dist(0.0,max_weight);
uniform_int_distribution<int> distint(0,num_particles-1);

int index = distint(gen);
double beta = 0.0;

// Implementing the wheel concept
vector<Particle> resampled;
for(int i=0;i<num_particles;i++)
{
 beta = beta + dist(gen)*2.0f;
 while(beta > weights[index])
 {
  beta= beta - weights[index];
  index = (index+1)%num_particles;
 } 
 resampled.push_back(particles[index]);
}

particles = resampled;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
