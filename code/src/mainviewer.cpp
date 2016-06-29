#include <mecacell/mecacell.h>
#include <mecacell/viewer/viewer.h>
#include "core/typesconfig.hpp"
#include "external/cxxopts.hpp"
#include "viewer/viewerextensions.hpp"

using G = GRN<RealCoords>;
using MG = MGRN<MGClassic>;
using C = GRNPlantController<G>;
using MC = GRNPlantController<MG>;
template <typename Ctrl> using PC = PlantCell<Ctrl, MecaCell::VolumeMembrane>;

template <typename V> int launchViewer(V& v) {
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

int main(int argc, char** argv) {
	std::string evaluatorName;
	std::string grnType = "grn";
	if (grnType == "grn") {
		MecacellViewer::Viewer<Scenario<PC<C>>> v(argc, argv);
		launchViewer(v);
	} else if (grnType == "mgrn") {
		MecacellViewer::Viewer<Scenario<PC<MC>>> v(argc, argv);
		launchViewer(v);
	} else {
		std::cerr << "No valid grn type found, aborting." << std::endl;
		return 1;
	}
	return 0;
}
