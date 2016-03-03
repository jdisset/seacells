#include <mecacell/mecacell.h>
#include <mecacell/viewer/viewer.h>
#include "core/typesconfig.hpp"
#include "viewer/grnviewer.hpp"
#include "viewer/viewerextensions.hpp"

int main(int argc, char *argv[]) {
	MecacellViewer::Viewer<TypesConfig::ScenarioType> v(argc, argv);
	MorphologyCapture mc;
	PewPew pp;
	ActiveConnectionsView ac;
	GroundView g;
	NutrientsView nv;
	GRNViewer<TypesConfig::GrnType> gv;
	v.registerPlugin(mc);
	v.registerPlugin(pp);
	v.registerPlugin(gv);
	v.registerPlugin(nv);
	v.registerPlugin(ac);
	v.registerPlugin(g);
	return v.exec();
}
