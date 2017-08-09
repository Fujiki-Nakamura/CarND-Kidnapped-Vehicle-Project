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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    default_random_engine gen;
    // create normal distribution for x, y and theta
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    // initialize particles
    num_particles = 100;
    particles.resize(num_particles);
    weights.resize(num_particles);
    double initial_weight = 1.0 / num_particles;

    for (int i = 0; i < num_particles; i++) {
        Particle p;
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        p.weight = initial_weight;
        paticles.push_back(p);
        weights.push_back(initial_weights);
    }
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    double velocity_delta_t = velocity * delta_t;
    double velocity_yaw_rate = velocity / yaw_rate;
    double yaw_rate_delta_t = yaw_rate * delta_t;

    default_random_engine gen;
    normal_distribution<double> dist_x(0.0, std_pos[0]);
    normal_distribution<double> dist_y(0.0, std_pos[1]);
    normal_distribution<double> dist_theta(0.0, std_pos[2]);

    for (int i = 0; i < num_particles; i++) {
        if (fabs(yaw_rate) < 0.001) {
            particles[i].x = velocity_delta_t * cos(paticles[i].theta);
            particles[i].y = velocity_delta_t * sin(paticles[i].theta);
        }
        else {
            particles[i].x = velocity_yaw_rate * (sin(particles[i].theta + yaw_rate_delta_t) - sin(paticles[i].theta));
            particles[i].y = velocity_yaw_rate * (cos(particles[i].theta) - cos(paticles[i].theta + yaw_rate_delta_t))
            particles[i].theta = particles[i].theta + yaw_rate_delta_t;
        }
        // add noises
        particles[i].x      += dist_x(gen);
        particles[i].y      += dist_y(gen);
        particles[i].theta  += dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    for (int i = 0; i < observations.size(); i++) {
        double min_distance = numeric_limits<double>::max();
        
        for (int j = 0; j < predicted.size; j++) {
            double distance = dist(
                observations[i].x, observations[i].y,
                predicted[j].x, predicted[j].y
            );
            if (distance < min_distance) {
                min_distance = distance;
                observations[i].id = predicted[j].id;
            }
        }
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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
    for (int i = 0; i < num_particles; i++) {
        double p_x      = particles[i].x;
        double p_y      = particles[i].y;
        double p_theta  = particles[i].theta;

        vector<LandmarkObs> predictions;
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            int landmark_id     = map_landmarks.landmark_list[j].id_i;
            double landmark_x   = map_landmarks.landmark_list[j].x_f;
            double landmark_y   = map_landmarks.landmark_list[j].y_f;

            double distance = dist(p_x, p_y, landmark_x, landmark_y);
            if (distance < sensor_range) {
                LandmarkObs landmark_in_range = {
                    landmark_id, landmark_x, landmark_y
                };
                predictions.push_back(landmark_in_range);
            }
        }

        vector<LandmarkObs> observations_map;
        double cos_theta = cos(p_theta);
        double sin_theta = sin(p_theta);
        for (int j = 0; j < observations.size(); j++) {
            LandmarkObs landmark_tmp;
            landmark_tmp.x = p_x + cos_theta * observations[j].x - sin_theta * observations[j].y;
            landmark_tmp.y = p_y + sin_theta * observations[j].x + cos_theta * observations[j].y;
            observations_map.push_back(landmark_tmp);
        }

        dataAssociation(predictions, observations_map);

        particles[i].weight = 1.0;
        for (int j = 0; j < observations_map.size(); j++) {
            int index = observations_map[j].id;
            double o_x = observations_map[j].x;
            double o_y = observations_map[j].y;
            
            double x_predicted = predictions[j].x;
            double y_predicted = predictions[j].y;

            double x_term = pow(o_x - x_predicted, 2) / (2 * pow(landmark[0], 2));
            double y_term = pow(o_y - y_predicted, 2) / (2 * pow(landmark[1], 2));
            double w = exp(-(x_term + y_term)) / (2 * M_PI * std_landmark[0] * std_landmark[1]);
            particles.weight *= w;
        }
        weights[i] = particles.weight;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    vector<Particle> particles_resampled;
    particles_resampled.resize(num_particles);
    default_random_engine gen;
    discrete_distribution<int> dist(weights.begin(), weights.end());
    
    for (int i = 0; i < num_particles; i ++) {
        int index = dist(gen);
        particles_resampled[i] = particles[index];
    }
    particles = particles_resampled;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
