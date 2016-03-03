#include <mecacell/mecacell.h>
#include "core/typesconfig.hpp"

int main(int argc, char *argv[]) {
	TypesConfig::ScenarioType sc;
	sc.init(argc, argv);
	while (!sc.finished()) {
		sc.loop();
		std::cout << "updt " << sc.getWorld().getNbUpdates() << ", "
		          << sc.getWorld().cells.size() << " cells" << std::endl;
		int i = 0;
		// for (auto &n : sc.nutrientSources) {
		// std::cerr << " N " << i++ << " : content = " << n.content << std::endl;
		//}
	}
	sc.terminate();
	return 0;
}
