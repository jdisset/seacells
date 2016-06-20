#ifndef TYPESCONFIG_HPP
#define TYPESCONFIG_HPP
#include <mecacell/mecacell.h>
#include "../external/grgen/common.h"
#include "../external/grgen/grn.hpp"
#include "../external/grgen/mgrn.hpp"
#include "../external/grgen/mgclassic.hpp"
#include "../external/grgen/real.hpp"
#include "config.hpp"
#include "plantcell.hpp"
#include "plantcontroller.hpp"
#include "scenario.hpp"

struct TypesConfig {
	// using GrnType = GRN<Classic>;
	// using CtrlType = GRNPlantController<GrnType>;
	// using CellType = PlantCell<CtrlType, MecaCell::VolumeMembrane>;
	//using ScenarioType = Scenario<CellType>;
};
#endif
