#include <integrator/integrator.hpp>
#include <iostream>
#include <chrono>

int main()
{
	int n = 3;
	int m = 3;
	int m0 = 4;
	int n0 = 3;
	int n1 = n0 + n + 1;
	int m1 = m0 + m + 1;
	int nm = (n + 1) * (m + 1);
    std::vector<std::vector<Fraction>> G;
    std::vector<std::vector<Fraction>>* CCC = nullptr;
		//G = []
		//relu = lambda x : max(x, 0)

	std::vector<std::vector<std::pair<Fraction, Fraction>>> traingles = 
    {
        { {Fraction(0), Fraction(0)}, {Fraction(0), Fraction(1,4)}, {Fraction(1,4), Fraction(1,4)} },
        { {Fraction(0), Fraction(0)}, {Fraction(1,4), Fraction(1,4)}, {Fraction(1,4), Fraction(0)} },
        { {Fraction(1,4), Fraction(0)}, {Fraction(1,4), Fraction(1,4)}, {Fraction(1,2), Fraction(1,4)} },
        { {Fraction(1,4), Fraction(0)}, {Fraction(1,2), Fraction(1,4)}, {Fraction(1,2), Fraction(0)} },
        { {Fraction(1,2), Fraction(0)}, {Fraction(1,2), Fraction(1,4)}, {Fraction(3,4), Fraction(1,4)} },
        { {Fraction(1,2), Fraction(0)}, {Fraction(3,4), Fraction(1,4)}, {Fraction(3,4), Fraction(0)} },
        { {Fraction(3,4), Fraction(0)}, {Fraction(3,4), Fraction(1,4)}, {Fraction(1), Fraction(1,4)} },
        { {Fraction(3,4), Fraction(0)}, {Fraction(1), Fraction(1,4)}, {Fraction(1), Fraction(0)} },
        { {Fraction(0), Fraction(1,4)}, {Fraction(0), Fraction(1,2)}, {Fraction(1,4), Fraction(1,2)} },
        { {Fraction(0), Fraction(1,4)}, {Fraction(1,4), Fraction(1,2)}, {Fraction(1,4), Fraction(1,4)} },
        { {Fraction(1,4), Fraction(1,4)}, {Fraction(1,4), Fraction(1,2)}, {Fraction(1,2), Fraction(1,2)} },
        { {Fraction(1,4), Fraction(1,4)}, {Fraction(1,2), Fraction(1,2)}, {Fraction(1,2), Fraction(1,4)} },
        { {Fraction(1,2), Fraction(1,4)}, {Fraction(1,2), Fraction(1,2)}, {Fraction(3,4), Fraction(1,2)} },
        { {Fraction(1,2), Fraction(1,4)}, {Fraction(3,4), Fraction(1,2)}, {Fraction(3,4), Fraction(1,4)} },
        { {Fraction(3,4), Fraction(1,4)}, {Fraction(3,4), Fraction(1,2)}, {Fraction(1), Fraction(1,2)} },
        { {Fraction(3,4), Fraction(1,4)}, {Fraction(1), Fraction(1,2)}, {Fraction(1), Fraction(1,4)} },
        { {Fraction(0), Fraction(1,2)}, {Fraction(0), Fraction(3,4)}, {Fraction(1,4), Fraction(3,4)} },
        { {Fraction(0), Fraction(1,2)}, {Fraction(1,4), Fraction(3,4)}, {Fraction(1,4), Fraction(1,2)} },
        { {Fraction(1,4), Fraction(1,2)}, {Fraction(1,4), Fraction(3,4)}, {Fraction(1,2), Fraction(3,4)} },
        { {Fraction(1,4), Fraction(1,2)}, {Fraction(1,2), Fraction(3,4)}, {Fraction(1,2), Fraction(1,2)} },
        { {Fraction(1,2), Fraction(1,2)}, {Fraction(1,2), Fraction(3,4)}, {Fraction(3,4), Fraction(3,4)} },
        { {Fraction(1,2), Fraction(1,2)}, {Fraction(3,4), Fraction(3,4)}, {Fraction(3,4), Fraction(1,2)} },
        { {Fraction(3,4), Fraction(1,2)}, {Fraction(3,4), Fraction(3,4)}, {Fraction(1), Fraction(3,4)} },
        { {Fraction(3,4), Fraction(1,2)}, {Fraction(1), Fraction(3,4)}, {Fraction(1), Fraction(1,2)} },
        { {Fraction(0), Fraction(3,4)}, {Fraction(0), Fraction(1)}, {Fraction(1,4), Fraction(3,4)} },
        { {Fraction(0), Fraction(3,4)}, {Fraction(1,4), Fraction(1)}, {Fraction(1,4), Fraction(1,2)} }, 
        { {Fraction(1,4), Fraction(3,4)}, {Fraction(1,4), Fraction(1)}, {Fraction(1,2), Fraction(1)} },
        { {Fraction(1,4), Fraction(3,4)}, {Fraction(1,2), Fraction(1)}, {Fraction(1,2), Fraction(3,4)} },
        { {Fraction(1,2), Fraction(3,4)}, {Fraction(1,2), Fraction(1)}, {Fraction(3,4), Fraction(1)} },
        { {Fraction(1,2), Fraction(3,4)}, {Fraction(3,4), Fraction(1)}, {Fraction(3,4), Fraction(3,4)} },
        { {Fraction(3,4), Fraction(3,4)}, {Fraction(3,4), Fraction(1)}, {Fraction(1), Fraction(1)} },
        { {Fraction(3,4), Fraction(3,4)}, {Fraction(1), Fraction(1)}, {Fraction(1), Fraction(3,4)} }
    };

    for (auto& triangle : traingles)
    {
        std::vector<std::vector<Fraction>> G_k;
        for (int i = 0; i < nm; i++)
        {
            std::vector<Fraction> temp;
            for (int j = 0; j < nm; j++)
            {
                temp.push_back(Fraction(0));
            }
            G_k.push_back(temp);
        }
        Integrator integrator = Integrator(triangle);
        if (CCC == nullptr)
        {
            CCC = integrator.Binome;
        }
        else
        {
            integrator.Binome = CCC;
        }
        auto start_time = std::chrono::high_resolution_clock::now();
        for (int i1 = n0; i1 < n1; i1++)
        {
            for (int j1 = m0; j1 < m1; j1++)
            {
                for (int i2 = n0; i2 < n1; i2++)
                {
                    for (int j2 = m0; j2 < m1; j2++)
                    {
                        Fraction S, S1, S2, S3;
                        if (i2 >= 4)
                        {
                            S1 = integrator.Integrate(i1 + i2 - 4, j1 + j2) * (i2 - 3) * (i2 - 2) * (i2 - 1) * (i2);
                        }
                        if (i2 >= 2 && j2 >= 2)
                        {
                            S2 = Fraction(2) * integrator.Integrate(i1 + i2 - 2, j1 + j2 - 2) * (i2) * (i2 - 1) * (j2 - 1) * (j2);
                        }
                        if (j2 >= 4)
                        {
                            S3 = integrator.Integrate(i1 + i2, j1 + j2 - 4) * (j2 - 3) * (j2 - 2) * (j2 - 1) * (j2);
                        }
                        S = S1 + S2 + S3;
                        G_k[(i1 - n0) * m + (j1 - m0)][(i2 - n0) * m + (j2 - m0)] = S;
                    }
                }
            }
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> ms_double = end_time - start_time;
        std::cout << "Iteration time: " << ms_double.count() << "ms/n";
        for (auto& vector : G_k)
        {
            for (auto& scalar : vector)
            {
                std::cout << scalar << ' ';
            }
            std::cout << '\n';
        }
    }

	return 0;
}