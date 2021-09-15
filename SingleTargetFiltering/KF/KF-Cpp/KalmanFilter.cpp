#include "KalmanFilter.hpp"

#include <cstring>  // header for the memcpy function:

#include <boost/numeric/ublas/matrix.hpp>

using boost::numeric::ublas::identity_matrix;
using boost::numeric::ublas::matrix;


const float INIT_Pk = 0.001;
const float KF_INITIAL_VELOCITY_x = 0; 
const float KF_INITIAL_VELOCITY_y = 0; 

const float sigma_v = 0.01;  // std of Q_k
const float sigma_r = 0.1;   // std of H_k

class KalmanFilter::Impl
{
public:
	explicit Impl(const float status[KF_SIZE]);
	const float* get_status();
	const float* get_prediction();
	bool update(const float status[KF_SIZE]);
	bool predict();
	bool skip();

private:
	matrix<float> Q_k = matrix<float>(KF_SIZE * 2, KF_SIZE * 2);  // process noise covariance
	matrix<float> H_k = matrix<float>(KF_SIZE, KF_SIZE * 2);      // observation model
	matrix<float> R_k = matrix<float>(KF_SIZE, KF_SIZE);          // observation noise covariance
	matrix<float> F_k = matrix<float>(KF_SIZE * 2, KF_SIZE * 2);  // state transition model

	matrix<float> x_k = matrix<float>(KF_SIZE * 2, 1);
	matrix<float> x_k_pred = matrix<float>(KF_SIZE * 2, 1);
	matrix<float> P_k = matrix<float>(KF_SIZE * 2, KF_SIZE * 2);
	matrix<float> P_k_pred = matrix<float>(KF_SIZE * 2, KF_SIZE * 2);
	
	bool order2_inverse(const matrix<float>& in, matrix<float>& out)
	{
		if (in.size1() != KF_SIZE || in.size2() != KF_SIZE)
		{
			return false;
		}

		float det = in(0, 0) * in(1, 1) - in(0, 1) * in(1, 0);
		if (det == 0)
		{
			return false;
		}

		out(0, 0) = in(1, 1) / det;
		out(1, 1) = in(0, 0) / det;
		out(0, 1) = in(0, 1) / det * (-1);
		out(1, 0) = in(1, 0) / det * (-1);

		return true;
	}

	float m_status[KF_SIZE * 2];
};

KalmanFilter::Impl::Impl(const float status[KF_SIZE])
{
	memcpy(m_status, status, sizeof(float) * KF_SIZE);
	for (size_t i = KF_SIZE; i < KF_SIZE * 2; ++i)
	{
		m_status[i] = 0.0f;
	}

	identity_matrix<float> tmp(KF_SIZE * 2);
	identity_matrix<float> tmp_r(KF_SIZE);

	/*  F_k  =  1 0 1 0
	 *          0 1 0 1
	 *          0 0 1 0
	 *          0 0 0 1
	 */
	F_k = tmp;
	F_k(0, 2) = 1;
	F_k(1, 3) = 1;

	/*  Q_k  =  1 0 0 0
	 *          0 1 0 0
	 *          0 0 1 0   *  sigma_v^2
	 *          0 0 0 1
	 */
	
	Q_k = pow(sigma_v,2) * tmp;

	/*   H_k = 1 0 0 0
	 *         0 1 0 0
	 */
	H_k.clear();
	H_k(0, 0) = 1;
	H_k(1, 1) = 1;

	/*   R_k = 1 0   * *  sigma_r^2
	 *         0 1
	 */
	R_k.clear();
	R_k = pow(sigma_r,2) * tmp_r;

	/*   H_k = 1 0 0 0
	 *         0 1 0 0
	 */
	H_k.clear();
	H_k(0, 0) = 1;
	H_k(1, 1) = 1;

	/*  x_k = (x,y,5,5)'
	 */
	x_k(0, 0) = status[0];
	x_k(1, 0) = status[1];
	x_k(2, 0) = KF_INITIAL_VELOCITY_x;
	x_k(3, 0) = KF_INITIAL_VELOCITY_y;


	/* P_k = I * INIT_Pk
	 */
	P_k = tmp * INIT_Pk;


	x_k_pred = x_k;
	P_k_pred = P_k;
}

const float* KalmanFilter::Impl::get_status()
{
	for (size_t i = 0; i < KF_SIZE * 2; ++i)
	{
		m_status[i] = x_k(i, 0);
	}
	return m_status;
}

const float* KalmanFilter::Impl::get_prediction()
{
	for (size_t i = 0; i < KF_SIZE * 2; ++i)
	{
		m_status[i] = x_k_pred(i, 0);
	}
	return m_status;
}

bool KalmanFilter::Impl::update(const float status[KF_SIZE])
{
	matrix<float> y_k(KF_SIZE, 1);  // Build y_k as same size of z_p
	y_k(0, 0) = status[0];
	y_k(1, 0) = status[1];
	y_k -= prod(H_k, x_k_pred);

	matrix<float> t1 = prod(H_k, P_k_pred);
	matrix<float> t2 = prod(t1, trans(H_k));
	matrix<float> S_k = t2 + R_k;

	matrix<float> K_k(P_k_pred);
	K_k = prod(K_k, trans(H_k));

	matrix<float> invS(S_k);
	if (!order2_inverse(S_k, invS))
	{
		return false;
	}

	K_k = prod(K_k, invS);

	// X(k+1)  = X(k-pred)+K(k)Y(k)
	x_k = x_k_pred + prod(K_k, y_k);

	// P(k+1)  = (I-K(k)H(k))P(k-pred)
	P_k = prod(K_k, H_k);
	identity_matrix<float> I(P_k.size1());
	P_k = I - P_k;
	P_k = prod(P_k, P_k_pred);
	return true;
}

bool KalmanFilter::Impl::predict()
{
	matrix<float> x_prev(x_k);
	matrix<float> P_prev(P_k);

	x_k_pred = prod(F_k, x_prev);
	P_k_pred = prod(F_k, P_prev);
	P_k_pred = prod(P_k_pred, trans(F_k)) + Q_k;
	return true;
}

bool KalmanFilter::Impl::skip()
{
	x_k = x_k_pred;
	P_k = P_k_pred;
	return true;
}

KalmanFilter::KalmanFilter(const float status[KF_SIZE]) : m_impl(new Impl(status))
{
}

KalmanFilter::~KalmanFilter() = default;

const float* KalmanFilter::get_status()
{
	return m_impl->get_status();
}

const float* KalmanFilter::get_prediction()
{
	return m_impl->get_prediction();
}

bool KalmanFilter::update(const float status[KF_SIZE])
{
	return m_impl->update(status);
}
bool KalmanFilter::predict()
{
	return m_impl->predict();
}

bool KalmanFilter::skip()
{
	return m_impl->skip();
}

