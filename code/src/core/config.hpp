#ifndef CONFIG_HPP
#define CONFIG_HPP
#include <array>
#define WATER 0
#define LIGHT 1

struct Config {
	// init
	static constexpr double EPSILON_GROUND = 50.0;
	static constexpr double STEMCELL_Y = -60.0;
	static constexpr double STEMCELL_NUT0 = 0.15;
	static constexpr double STEMCELL_NUT1 = 0.5;

	// world caracs
	static constexpr double SIM_DT = 1.0 / 150.0;
	static constexpr unsigned int GRN_STEPS_PER_UPDATE = 1;
	static constexpr double AIR_VISCOSITY = 0.0005;
	static constexpr double GROUND_VISCOSITY = 0.001;
	static constexpr double GRAVITY = -10.0;
	static constexpr double GROUND_REACTION = 0.0002;

	// evo
	static constexpr unsigned int INITIAL_NB_REGULS = 1;
	static constexpr double DEFAULT_SIM_DURATION = 500.0;
	static constexpr double DEFAULT_MAX_CELLS = 420;

	// cells caracs
	static constexpr double CELL_GROWTH_SPEED = 0.35;
	static constexpr bool ENABLE_SOLIDIFY = true;
	static constexpr double MIN_NUTRIENTS_FOR_DIVISION = 0.0;
	static constexpr double DIVISION_NRJ_CONSUMPTION = 0.015;
	static constexpr double NORMAL_NRJ_CONSUMPTION = 0.012;
	static constexpr unsigned int NB_MORPHOGENS = 3;
	static const std::array<double, 5> morphoDiffusionCoefs;
	static constexpr double MORPHOGEN_UPDATE_INTERVAL = 5.0 * SIM_DT;
	static constexpr double MORPHO_SAMPLING_DIST = 40.0;
	static constexpr double NUTRIENT_SAMPLING_DIST = 40.0;
	static constexpr double NUTRIENT_SAMPLING_COEF = 3.0;

	// Nutrients
	static constexpr unsigned int NB_NUTRIENTS = 2;
	static constexpr unsigned int NB_NUTRIENTS_SOURCES = 185;
	static constexpr double NUTRIENTS_BOUNDING_AREA = 1700.0;
	static constexpr double NUTRIENTS_BOUNDING_DEPTH = 1000.0;
	static constexpr double TYPICAL_NUTRIENTS_RADIUS = 240.0;
	static constexpr double NUTRIENTS_VISCOSITY = 1.0;
	static constexpr double NUTRIENTS_DIFFUSION_K = 0.045;
	static constexpr double NUTRIENT_QUANTITY = 0.035;
	static constexpr double NUTRIENT_DEPTH_INCREASE_COEF = 0.033;
	static constexpr double NUTRIENT_DEPTH_INCREASE_POW = 1.45;
	static constexpr double FIRST_NUTRIENT_SOURCE_THRESHOLD = 100.0;
	static constexpr double SUN_INTENSITY = 0.5;
	static constexpr double MAX_LIGHT_THRESHOLD = 300.0;
};
#endif
