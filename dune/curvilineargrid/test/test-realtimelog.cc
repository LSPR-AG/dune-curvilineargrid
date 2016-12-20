#include <config.h>

#include <cstdlib>
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/common/realtimelog.hh>

using namespace Dune;

using namespace Dune::CurvGrid;


int main (int argc , char **argv) {
	static MPIHelper & mpihelper = MPIHelper::instance(argc, argv);

	CurvGridRealTimeLog::init(mpihelper, "memlog", 1, 5);


	CurvGridRealTimeLog::note("Started Math");

	/* initialize random seed: */
	srand (time(NULL) * mpihelper.rank());
	{
		const int WORK_SIZE = 1000000000;
		const int LOG_SIZE = 100;

		// Do some whatever n^2 algorithm that takes infinite amount of memory
		std::vector<double> m {1.0/3, 1.0/7};
		int work = floor(WORK_SIZE * (double(rand()) / RAND_MAX));

		std::cout << "Suggested work " << work << std::endl;

		for (int i = 0; i < work; i++) {
			m.push_back(i * 0.5);

			if (i % (WORK_SIZE / LOG_SIZE) == 0 ) {
				CurvGridRealTimeLog::logvar("workvar", std::to_string(m[i]));
			}
		}
	}

	CurvGridRealTimeLog::note("Finished Math");

	std::this_thread::sleep_for(std::chrono::seconds(5));

	CurvGridRealTimeLog::stop();

	return 0;
}
