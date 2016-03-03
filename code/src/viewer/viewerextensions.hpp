#ifndef VIEWEREXTENSIONS_HPP
#define VIEWEREXTENSIONS_HPP
#include <QString>
#include <QOpenGLShaderProgram>
#include <iostream>
#include <random>
#include "drawextensions.hpp"
#include "../core/capture.hpp"

using namespace MecacellViewer;

struct ForcesControls {
 public:
	template <typename R> void onLoad(R* renderer) {
		renderer->addButton("rdmForce", "GENERALACTIONS_MENU", "Random force",
		                    [=](R* r, Button<R>*) {
			                    for (auto& c : r->getScenario().getWorld().cells) {
				                    c->receiveExternalForce(MecaCell::Vec::randomUnit() * 1000.0);
			                    }
			                  })->setColor(59, 190, 220, 150);
		renderer->addButton("rdmVelocity", "GENERALACTIONS_MENU", "Random Velocity",
		                    [=](R* r, Button<R>*) {
			                    for (auto& c : r->getScenario().getWorld().cells) {
				                    c->setVelocity(MecaCell::Vec::randomUnit() * 80.0);
			                    }
			                  })->setColor(59, 190, 200, 180);
		renderer->addButton("RandomTorque", "GENERALACTIONS_MENU", "Add torque",
		                    [=](R* r, Button<R>*) {
			                    if (r->getSelectedCell()) {
				                    r->getSelectedCell()->receiveExternalTorque(
				                        MecaCell::Vec(0, 1, 0) * 100000.0);
			                    }
			                  })->setColor(250, 110, 20, 150);
		renderer->addButton("RandomAngVel", "GENERALACTIONS_MENU", "Random Angular Velocity",
		                    [=](R* r, Button<R>*) {
			                    for (auto& c : r->getScenario().getWorld().cells) {
				                    c->setAngularVelocity(MecaCell::Vec::randomUnit() * 0.5);
			                    }
			                  })->setColor(230, 90, 50, 150);
		renderer->addButton(
		    "Gotocenter", "GENERALACTIONS_MENU", "Force toward (0,0,0)",
		    [=](R* r, Button<R>*) {
			    for (auto& c : r->getScenario().getWorld().cells) {
				    c->receiveExternalForce(-c->getPosition().deltaDirection(0.3).normalized() *
				                            4000.0);
			    }
			  });
	}
};

class PewPew {
 public:
	template <typename R> void onLoad(R* renderer) {
		using C = typename R::Cell;
		using V = typename R::Vec;
		renderer->addButton("addcell", "GENERALACTIONS_MENU", "Add cell",
		                    [=](R* r, Button<R>*) {
			                    V camPos(r->getCamera().getPosition());
			                    V camDir(r->getCamera().getOrientation());
			                    V camUp(r->getCamera().getUpVector());
			                    C* nc = new C(camPos + camDir * 50.0);
			                    r->getScenario().getWorld().addCell(nc);
			                  });
		renderer->addButton("pew", "GENERALACTIONS_MENU", "Pew!", [=](R* r, Button<R>*) {
			V camPos(r->getCamera().getPosition());
			V camDir(r->getCamera().getOrientation());
			V camUp(r->getCamera().getUpVector());
			C* nc = new C(camPos + camDir * 30.0 - camUp * 30.0);
			nc->setVelocity(camDir * 1500.0);
			r->getScenario().getWorld().addCell(nc);
		});
		Button<R>* b = renderer->addButton(
		    "pewpew", "GENERALACTIONS_MENU", "Pew Pew!", [=](R* r, Button<R>*) {
			    std::random_device rd;
			    std::mt19937 gen(rd());
			    std::uniform_real_distribution<double> dis(0.0, 1.0);
			    V camPos(r->getCamera().getPosition());
			    V camDir(r->getCamera().getOrientation());
			    V camUp(r->getCamera().getUpVector());
			    for (int i = 0; i < 30; ++i) {
				    V rnd = V::randomUnit() * dis(gen) * 100.0;
				    C* nc = new C(camPos + camDir * 30.0 - camUp * 30.0 + rnd);
				    nc->setVelocity(camDir * 2000.0);
				    r->getScenario().getWorld().addCell(nc);
			    }
			  });
		b->setColor(200, 20, 30, 100);
	}
};
class Division {
 public:
	template <typename R> void onLoad(R* renderer) {
		renderer->addButton("divide", "GENERALACTIONS_MENU", "Start division",
		                    [=](R* r, Button<R>*) {
			                    if (r->getSelectedCell()) {
				                    r->getSelectedCell()->startDivision();
			                    }
			                  })->setColor(123, 0, 250, 150);
	}
};

class MorphologyCapture {
 public:
	template <typename R> void onLoad(R* renderer) {
		renderer->addButton("morphcapt", "GENERALACTIONS_MENU", "Morphology capture",
		                    [=](R* r, Button<R>*) {
			                    auto capture = ClusterTools::getMatrixCapture(
			                        r->getScenario().getWorld().cells);
			                    std::cerr << " ---- ---- ---- ---- " << std::endl;
			                    std::cerr << ClusterTools::captMatrixToString(capture, true);
			                  })->setColor(123, 120, 250, 150);
	}
};

class DetailedConnectionsView {
	DetailedConnections ac;

 public:
	template <typename R> void onLoad(R* renderer) {
		MenuElement<R>* nativeDisplayMenu = renderer->getDisplayMenu();
		MenuElement<R> contactView = {"Detailed connections", false};
		ac.load();
		contactView.onToggled = [&](R* r, MenuElement<R>* me) {
			if (me->isChecked()) {
				r->addPaintStepsMethods(19, [&](R* r2) { ac.call(r2); });
			} else {
				r->erasePaintStepsMethods(19);
			}
		};
		nativeDisplayMenu->at("Cells").add(contactView);
	}
};

class ActiveConnectionsView {
	ActiveConnections ac;

 public:
	template <typename R> void onLoad(R* renderer) {
		MenuElement<R>* nativeDisplayMenu = renderer->getDisplayMenu();
		MenuElement<R> contactView = {"Active connections", false};
		ac.load();
		contactView.onToggled = [&](R* r, MenuElement<R>* me) {
			if (me->isChecked()) {
				r->addPaintStepsMethods(18, [&](R* r2) { ac.call(r2); });
			} else {
				r->erasePaintStepsMethods(18);
			}
		};
		nativeDisplayMenu->at("Cells").add(contactView);
	}
};

class GroundView {
	Ground ac;

 public:
	template <typename R> void onLoad(R* renderer) {
		MenuElement<R>* nativeDisplayMenu = renderer->getDisplayMenu();
		MenuElement<R> contactView = {"Ground", false};
		ac.load();
		contactView.onToggled = [&](R* r, MenuElement<R>* me) {
			if (me->isChecked()) {
				r->addPaintStepsMethods(50, [&](R* r2) { ac.call(r2); });
			} else {
				r->erasePaintStepsMethods(50);
			}
		};
		nativeDisplayMenu->at("Cells").add(contactView);
	}
};

class NutrientsView {
	Nutrients nv;

 public:
	template <typename R> void onLoad(R* renderer) {
		MenuElement<R>* nativeDisplayMenu = renderer->getDisplayMenu();
		MenuElement<R> contactView = {"Nutrients", false};
		nv.load();
		contactView.onToggled = [&](R* r, MenuElement<R>* me) {
			if (me->isChecked())
				r->addPaintStepsMethods(40, [&](R* r2) { nv.call(r2); });
			else
				r->erasePaintStepsMethods(40);
		};
		nativeDisplayMenu->at("Cells").add(contactView);
	}
};

#endif
