#include <mecacell/mecacell.h>
#include <chrono>
#include "core/typesconfig.hpp"

using G = GRN<RealCoords>;
using MG = MGRN<MGClassic>;
using C = GRNPlantController<G>;
using MC = GRNPlantController<MG>;
template <typename Ctrl> using PC = PlantCell<Ctrl, MecaCell::VolumeMembrane>;

int main(int argc, char *argv[]) {
	Scenario<PC<C>> sc;
	auto start = std::chrono::system_clock::now();
	sc.init(argc, argv);
	while (!sc.finished()) {
		sc.loop();
	}
	sc.terminate();
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> diff = end - start;
	std::cerr << "Evaluated in " << diff.count() << "s.  Survival time = "
	          << sc.getWorld().getNbUpdates() * sc.getWorld().getDt() << std::endl;
	return 0;
}
