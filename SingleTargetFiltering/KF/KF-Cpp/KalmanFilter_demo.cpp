#include <iostream>
#include <math.h>
#include <cstring>  // header for the memcpy function:
#include<fstream>

#include "KalmanFilter.hpp"

using namespace std;


// int main (int argc, char **argv)   // the same
int main(int argc, char *argv[])
{
	float ob_noise = 20.0f;
	if (argc >= 2)
	{
		ob_noise = strtof(argv[1], nullptr);
		cout << "ob_noise = " << ob_noise << "\n";
	}

	srand(time(nullptr));

	float angle = 45 / M_PI / 2;
	float speed = 7.07f;
	float acc = 0.0f;
	float delta_angle = 1 / M_PI / 2;
	float observation[2] = {0, 0};

	KalmanFilter kf(observation);

	float ybx = 0;
	float yby = 0;
	const size_t test_len = 100;  // 50, 100
	float input[test_len][2];
	float output[test_len][4];

	for (size_t i = 0; i < test_len; ++i)
	{
		kf.predict();

		// Returns a 4x1 vector: x coord, y coord, x speed, y speed
		const float *xk = kf.get_status();
		memcpy(input[i], observation, sizeof(float) * 2);
		memcpy(output[i], xk, sizeof(float) * 4);

		observation[0] += ybx + static_cast<float>(rand()) / RAND_MAX * ob_noise * 2 - ob_noise;
		observation[1] += yby + static_cast<float>(rand()) / RAND_MAX * ob_noise * 2 - ob_noise;

		kf.update(observation);

		speed += acc;
		angle += delta_angle;
		ybx += speed * std::cos(angle);
		yby += speed * std::sin(angle);
	}

	cout << "\tInput\t\tOutput\n";
	cout.precision(2);
	cout << fixed;  // std::fixed

	for (size_t i = 0; i < test_len; ++i)
	{
		cout << i << ":\t";

		for (float in : input[i])
		{
			cout << in << ", ";
		}

		cout << "\t";

		for (float out : output[i])
		{
			cout << out << ", ";
		}

		cout << "\n";
	}

	// Save the input and the output
	std::string OutputFilename("../output.txt");  // filename of output data
	std::ofstream Out(OutputFilename);
	for (size_t i = 0; i < test_len; ++i)
	{
		Out << input[i][0] << " " << input[i][1] << " " << output[i][0] << " " << output[i][1] << endl;  // Ignore speed for display
	}
	Out.close();


	cout << "Done \n";
}



