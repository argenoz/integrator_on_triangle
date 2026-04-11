#include <integrator/integrator.hpp>

Integrator::Integrator(std::vector<std::pair<Fraction, Fraction>> triangle)
{
	ax = Fraction(triangle[0].first);
	ay = Fraction(triangle[0].second);
	bx = Fraction(triangle[1].first);
	by = Fraction(triangle[1].second);
	cx = Fraction(triangle[2].first);
	cy = Fraction(triangle[2].second);

	Kx = bx - cx;
	Ky = by - cy;
	acx = ax - cx;
	acy = ay - cy;
	J = acx * Ky - acy * Ky;
	std::vector<Fraction> First;
	First.push_back(Fraction(1));
	std::vector<Fraction> Second;
	Second.push_back(Fraction(1));
	Second.push_back(Fraction(1));

	Binome = new std::vector<std::vector<Fraction>>;
	Binome->push_back(First);
	Binome->push_back(Second);
}

Fraction Integrator::Cnk(int n, int k)
{
	while (n >= Binome->size())
	{
		std::vector<Fraction> last_layer = Binome->back();
		std::vector<Fraction> new_layer;
		new_layer.push_back(1);
		for (size_t i = 1; i < last_layer.size(); i++)
		{
			new_layer.push_back(last_layer[i - 1] + last_layer[i]);
		}
		new_layer.push_back(1);
		Binome->push_back(new_layer);
	}
	return (*Binome)[n][k];
}

Fraction Integrator::Integrate(int n, int m)
{
	Fraction Integral_a;
	Fraction Integral_b;
	Fraction Integral_c;
	Fraction Integral_d;
	Fraction Integral_e;
	for (int p = 0; p <= n; p++)
	{
		Integral_b = Fraction(0);
		int n_p = n - p;
		Fraction temp_ = Kx;
		for (int temp = 0; temp < p; temp++)
		{
			temp_ = temp_ * temp_;
		}
		Fraction theta1 = temp_ * Cnk(n, p);
		for (int q = 0; q <= m; q++)
		{
			Integral_c = Fraction(0);
			temp_ = Ky;
			for (int temp = 0; temp < q; temp++)
			{
				temp_ = temp_ * temp_;
			}
			Fraction theta = theta1 * temp_ * Cnk(n, p);
			int m_q = m - q;
			Fraction pq1 = Fraction(p + q + 1);
			for (int v = 0; v <= n_p; v++)
			{
				Integral_d = Fraction(0);
				temp_ = acx;
				Fraction temp__ = cx;
				for (int temp = 0; temp < v; temp++)
				{
					temp_ = temp_ * temp_;
				}
				for (int temp = 0; temp < n_p - v; temp++)
				{
					temp__ = temp__ * temp__;
				}
				Fraction La1 = temp_ * temp__ * Cnk(n_p, v);
				for (int u = 0; u <= m_q; u++)
				{
					int uv1 = v + u + 1;
					Integral_e = Fraction(0);
					temp_ = acy;
					temp__ = cy;
					for (int temp = 0; temp < u; temp++)
					{
						temp_ = temp_ * temp_;
					}
					for (int temp = 0; temp < m_q - u; temp++)
					{
						temp__ = temp__ * temp__;
					}
					Fraction La = La1 * temp_ * temp__ * Cnk(m_q, u) / Fraction(uv1);
					for (int s = 0; s <= uv1; s++)
					{
						temp_ = Cnk(uv1, s) / (pq1 + s);
						if (s % 2 == 0)
						{
							Integral_e = Integral_e + temp_;
						}
						else
						{
							Integral_e = Integral_e - temp_;
						}
					}
					Integral_d = Integral_d + Integral_e * La;
				}
				Integral_c = Integral_c + Integral_d;
			}
			Integral_b = Integral_b + Integral_c * theta;
		}
		Integral_a = Integral_a + Integral_b;
	}
	return Integral_a * J;
}
