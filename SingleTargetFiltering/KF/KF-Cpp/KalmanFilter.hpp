#pragma once

#include <memory>


const int KF_SIZE = 2;

class KalmanFilter 
{
public:
	explicit KalmanFilter(const float status[2]);
	~KalmanFilter();

	// Returns a 4x1 vector: x coord, y coord, x speed, y speed
	const float* get_status();
	const float* get_prediction();
	bool update(const float status[2]);
	bool predict();
	bool skip();

private:
	class Impl;
	std::unique_ptr<Impl> m_impl;
};

