#ifndef TYPESCONFIG_HPP
#define TYPESCONFIG_HPP
#include "../external/grgen/common.h"
#include "../external/grgen/grn.hpp"
#include "../external/grgen/classic.hpp"
#include "plantcell.hpp"
#include "plantcontroller.hpp"
#include "scenario.hpp"
#include "config.hpp"
#include <mecacell/mecacell.h>

struct TypesConfig {
	using GrnType = GRN<Classic>;
	using CtrlType = GRNPlantController<GrnType>;
	using CellType = PlantCell<CtrlType, MecaCell::VolumeMembrane>;
	using ScenarioType = Scenario<CellType>;
};
#endif
