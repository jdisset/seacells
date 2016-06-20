#include <mecacell/mecacell.h>
#include <mecacell/viewer/viewer.h>
#include "core/typesconfig.hpp"
#include "external/cxxopts.hpp"
#include "viewer/grnviewer.hpp"
#include "viewer/viewerextensions.hpp"
using G = GRN<RealCoords>;
using MG = MGRN<MGClassic>;
using C = GRNPlantController<G>;
using MC = GRNPlantController<MG>;

template <typename Ctrl> int launchViewer(int argc, char** argv) {
	MecacellViewer::Viewer<Scenario<PlantCell<Ctrl, MecaCell::VolumeMembrane>>> v(argc,
	                                                                              argv);
	MorphologyCapture mc;
	PewPew pp;
	ActiveConnectionsView ac;
	GroundView g;
	NutrientsView nv;
	v.registerPlugin(mc);
	v.registerPlugin(pp);
	v.registerPlugin(nv);
	v.registerPlugin(ac);
	v.registerPlugin(g);
	return v.exec();
}

int main(int argc, char* argv[]) {
	std::string grnType;
	try {
		cxxopts::Options options(argv[0]);
		options.add_options()("g,grn", "grn type", cxxopts::value<std::string>(grnType));
		options.parse(argc, argv);
	} catch (const cxxopts::OptionException& e) {
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	} catch (const std::bad_cast& e) {
		std::cout << "bad cast: " << e.what() << std::endl;
		exit(1);
	}

	if (grnType == "grn") {
		return launchViewer<C>(argc, argv);
	} else if (grnType == "mgrn") {
		return launchViewer<MC>(argc, argv);
	}
	std::cerr << "No valid grn type found, aborting." << std::endl;
	return 0;
}
